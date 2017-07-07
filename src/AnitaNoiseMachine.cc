#include "AnitaNoiseMachine.h"


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaNoiseMachine::AnitaNoiseMachine() {
  zeroInternals();
}




void AnitaNoiseMachine::zeroInternals() {

  rmsFifoPos = 0;
  rmsFifoFillFlag = false;
  memset(rmsFifo,0,NUM_PHI*NUM_ANTENNA_RINGS*NUM_POLS*fifoLength*sizeof(double)); 

  memset(avgMaps,0,NUM_POLS*nPhi*nTheta*sizeof(double));

  for (int poli=0; poli<NUM_POLS; poli++) {
    for (int loc=0; loc<fifoLength; loc++) {
      if (mapFifo[poli][loc] != NULL) {
	delete mapFifo[poli][loc];
	mapFifo[poli][loc] = NULL;
      }
    }
  }
  mapFifoPos = 0;
  rmsFifoFillFlag = false;

}

/*=======================
  Waveform RMS stuff */
void AnitaNoiseMachine::fillAvgRMSNoise(FilteredAnitaEvent *filtered) {
  
  rmsFifoPos++;
  if (rmsFifoPos >= fifoLength) {
    rmsFifoPos = 0;
    rmsFifoFillFlag = true;
  }

  for (int phi=0; phi<NUM_PHI; phi++) {
    for (int ringi=0; (AnitaRing::AnitaRing_t)ringi != AnitaRing::kNotARing; ringi++) {
      AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
      for (int poli=0; (AnitaPol::AnitaPol_t)poli != AnitaPol::kNotAPol; poli++) {
	AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)poli;

	const TGraphAligned *currWave = filtered->getRawGraph(phi,ring,pol)->even();
	double value = currWave->GetRMS(1);
	rmsFifo[phi][ringi][poli][rmsFifoPos] = currWave->GetRMS(1);
      }
    }
  }
}


double AnitaNoiseMachine::getAvgRMSNoise(int phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol){

  double value2 = 0;

  int endPoint = fifoLength;

  for (int pos=0; pos<endPoint; pos++) {
    value2 += pow(rmsFifo[phi][(int)ring][(int)pol][pos],2)/endPoint;
  }
  return sqrt(value2);

}
/*---------------------*/
	

/*======================
  Get the correlation maps from the analyzer and update everything internally */
void AnitaNoiseMachine::fillAvgMapNoise(UCorrelator::Analyzer *analyzer) {
  
  mapFifoPos++; //move to the next point for writing  

  //if you haven't even written once though stay at zero
  if (mapFifo[0][mapFifoPos] == NULL) mapFifoPos = 0;

  //wrap if you're at the end of the fifo
  if (mapFifoPos >= fifoLength) {
    mapFifoPos = 0;
    mapFifoFillFlag = true;
  }
  
  for (int poli=0; poli<NUM_POLS; poli++) {
    //delete old histogram if it is still there
    if (mapFifo[poli][mapFifoPos] != NULL) {
      delete mapFifo[poli][mapFifoPos];
      mapFifo[poli][mapFifoPos] = NULL;
    }
    
    //syntactically weirdly copy it out of UCorrelator    
    const UCorrelator::gui::Map *currMap = analyzer->getCorrelationMap((AnitaPol::AnitaPol_t)poli);
    TH2D newMap;
    currMap->Copy(newMap);
    mapFifo[poli][mapFifoPos] = (TH2D*)newMap.Clone();
  }

  //do the simple array too
  updateAvgMapNoise();

  return;

}

void AnitaNoiseMachine::updateAvgMapNoise() {
  /* This might be a more efficient way to do it since building a new giant average TProfile every time you want
     to find out what the noise was at a point is pretty lousy */

  //stolen from Rene
  for (int poli=0; poli<NUM_POLS; poli++) {
    for (Int_t iPhi = 0; iPhi < nPhi; iPhi++) {
      for (Int_t iTheta = 0; iTheta < nTheta; iTheta++) {     
	//subtract fifo position that is expiring (if there is one)
	if (mapFifoFillFlag) {
	  if (mapFifoPos == 0) 
	    { avgMaps[poli][iPhi][iTheta] -= mapFifo[poli][fifoLength-1]->GetBinContent(iPhi,iTheta); }
	  else
	    { avgMaps[poli][iPhi][iTheta] -= mapFifo[poli][mapFifoPos-1]->GetBinContent(iPhi,iTheta); }
	}
	
	//and add the new addition to the fifo
	avgMaps[poli][iPhi][iTheta] += mapFifo[poli][mapFifoPos]->GetBinContent(iPhi,iTheta);
      
      }
    }
  }
  
  
  return;

}

/*------------*/

TProfile2D* AnitaNoiseMachine::getAvgMapNoiseProfile(AnitaPol::AnitaPol_t pol) {

  int nBinX = mapFifo[(int)pol][0]->GetNbinsX();
  int nBinY = mapFifo[(int)pol][0]->GetNbinsY();
  int xMin  = mapFifo[(int)pol][0]->GetXaxis()->GetBinLowEdge(1);
  int yMin  = mapFifo[(int)pol][0]->GetYaxis()->GetBinLowEdge(1);
  int xMax  = mapFifo[(int)pol][0]->GetXaxis()->GetBinUpEdge(nBinX);
  int yMax  = mapFifo[(int)pol][0]->GetYaxis()->GetBinUpEdge(nBinY);

  std::stringstream name;
  name.str("");
  name << "mapAvg_" << pol;
  
  TProfile2D *outProfile  = new TProfile2D(name.str().c_str(),"Average Interferometric Map", nBinX, xMin, xMax, nBinY, yMin, yMax);

  int lengthToDo = fifoLength;
  if (mapFifoFillFlag == false) lengthToDo = mapFifoPos;

  //stolen from Rene
  for (Int_t ix = 1; ix <= nBinX; ix++) {
    Double_t x = outProfile->GetXaxis()->GetBinCenter(ix);
    for (Int_t iy = 1; iy <= nBinY; iy++) {     
      Double_t y = outProfile->GetYaxis()->GetBinCenter(iy);
      for (Int_t i=0; i<lengthToDo; i++) {
	outProfile->Fill(x,y,mapFifo[(int)pol][i]->GetBinContent(ix,iy),1);
      }
    }
  }

  return outProfile;
}




/*===========================
  Copy the relevant things into the noise summary */
void AnitaNoiseMachine::fillNoiseSummary(AnitaNoiseSummary *noiseSummary) {
  if (!rmsFifoFillFlag) {
    std::cout << "WARNING in AnitaNoiseMachine::fillNoiseSummary(): Fifo hasn't been filled entirely yet, ";
    std::cout << " only gotten " << rmsFifoPos << " out of " << fifoLength << " values" << std::endl;
  }
  
  noiseSummary->fifoLength = fifoLength;
  noiseSummary->mapFifoFillFlag = mapFifoFillFlag;
  noiseSummary->rmsFifoFillFlag = rmsFifoFillFlag;


  if (fillAvgMap || fillMap) {
    for (int poli=0; (AnitaPol::AnitaPol_t)poli != AnitaPol::kNotAPol; poli++) {
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)poli;
      if (noiseSummary->avgMapProf[poli] != NULL) {
	delete noiseSummary->avgMapProf[poli];
	noiseSummary->avgMapProf[poli] = NULL;
      }
      if (fillMap) noiseSummary->avgMapProf[poli] = (TH2D*)mapFifo[poli][mapFifoPos]->Clone(); //non-average is default
      else noiseSummary->avgMapProf[poli] = getAvgMapNoiseProfile(pol);
    }
  }
 
  return;

}

/*-----------------------------*/



/*=========================
  Added a thing to AnitaEventSummary that this machine needs to fill */


void AnitaNoiseMachine::fillEventSummary(AnitaEventSummary *eventSummary) {

  if (mapFifo[0][0] == NULL) return; //can't do it if you haven't filled anything yet

  for (int poli=0; poli<NUM_POLS; poli++) {
    for (int dir=0; dir<AnitaEventSummary::maxDirectionsPerPol; dir++) {
      double peakPhi = eventSummary->peak[poli][dir].phi;
      double peakTheta = eventSummary->peak[poli][dir].theta;
      int binx = mapFifo[poli][0]->GetXaxis()->FindBin(peakPhi);
      int biny = mapFifo[poli][0]->GetYaxis()->FindBin(peakTheta);
      double avgNoise = avgMaps[poli][binx][biny];
      eventSummary->peak[poli][dir].peakPastRMS = avgNoise;
    }
  }

  return;
}
