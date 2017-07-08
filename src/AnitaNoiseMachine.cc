#include "AnitaNoiseMachine.h"


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaNoiseMachine::AnitaNoiseMachine(const int length)
   : fifoLength(length)
{
  
  if (fifoLength < 2) {
    std::cout << "Warning from constructor of AnitaNoiseMachine(int length): ";
    std::cout << "You picked an average less than 2, this is a bad idea, and values will likely be wrong" << std::endl;
  }

  rmsFifo = (double*)malloc(NUM_PHI*NUM_ANTENNA_RINGS*NUM_POLS*fifoLength*sizeof(double));

  rollingMapAvg = (double*)malloc(NUM_POLS*nPhi*nTheta*sizeof(double));

  for (int poli=0; poli<NUM_POLS; poli++) {
    mapFifo[poli] = (TH2D**)malloc(fifoLength*sizeof(TH2D*));
    for (int fifoPos=0; fifoPos<fifoLength; fifoPos++) {
      mapFifo[poli][fifoPos] = NULL;
    }
  }
  
  //maps need to be initialized to NULL or the reset wont work
  for (int poli=0; poli<NUM_POLS; poli++) {
    for (int fifoPos=0; fifoPos<fifoLength; fifoPos++) {
      mapFifo[poli][fifoPos] = NULL;
    }
  }

  zeroInternals();
}

//so that I always am refering to the same index
int AnitaNoiseMachine::rmsFifoIndex(int phi, int ringi, int poli, int fifoLength) {
  return phi*(NUM_ANTENNA_RINGS*NUM_POLS*fifoLength) + ringi*(NUM_POLS*fifoLength) + poli*(fifoLength) + fifoLength;
}

//so that I always am refering to the same index
int AnitaNoiseMachine::rollingMapIndex(int poli, int iPhi, int iTheta) {
  return poli*(nPhi*nTheta) + iPhi*(nTheta) + iTheta;
}


//Resetting to initial state
void AnitaNoiseMachine::zeroInternals() {

  fJustInitialized = true;

  //reset save flags
  fillArray = false;
  fillMap = false;
  fillAvgMap = false;

  //reset rms fifo
  rmsFifoPos = 0;
  rmsFifoFillFlag = false;
  memset(rmsFifo,0,NUM_PHI*NUM_ANTENNA_RINGS*NUM_POLS*fifoLength*sizeof(double)); 

  //reset map double array fifo
  memset(rollingMapAvg,0,NUM_POLS*nPhi*nTheta*sizeof(double));

  //reset map histogram fifo
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



/*===============
  Fill everything */

void AnitaNoiseMachine::updateMachine(UCorrelator::Analyzer *analyzer,FilteredAnitaEvent *filtered) {
  updateAvgRMSFifo(filtered);
  updateAvgMapFifo(analyzer);

  if (fJustInitialized) fJustInitialized = false;
  return;
}









/*=======================
  Waveform RMS stuff */
void AnitaNoiseMachine::updateAvgRMSFifo(FilteredAnitaEvent *filtered) {
  
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

	const TGraphAligned *currWave = filtered->getFilteredGraph(phi,ring,pol)->even();
	double value = currWave->GetRMS(1);
	rmsFifo[rmsFifoIndex(phi,ringi,poli,rmsFifoPos)] = currWave->GetRMS(1);
      }
    }
  }
}


double AnitaNoiseMachine::getAvgRMSNoise(int phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol){

  double value2 = 0;

  int endPoint = fifoLength;

  for (int pos=0; pos<endPoint; pos++) {
    value2 += pow(rmsFifo[rmsFifoIndex(phi,(int)ring,(int)pol,pos)],2)/endPoint;
  }
  return sqrt(value2);

}
/*---------------------*/
	

/*======================
  Get the correlation maps from the analyzer and update everything internally */
void AnitaNoiseMachine::updateAvgMapFifo(UCorrelator::Analyzer *analyzer) {
  
  mapFifoPos++; //move to the next point for writing    

  //if you haven't even written once though stay at zero
  if (fJustInitialized) mapFifoPos = 0;

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


  /* This might be a more efficient way to do it since building a new giant average TProfile every time you want
     to find out what the noise was at a point is pretty lousy */

  //stolen from Rene
  for (int poli=0; poli<NUM_POLS; poli++) {
    for (Int_t iPhi = 0; iPhi < nPhi; iPhi++) {
      for (Int_t iTheta = 0; iTheta < nTheta; iTheta++) {     

	//subtract fifo position that is expiring (if there is one)
	if (mapFifoFillFlag) {
	  int lastPos = mapFifoPos-1;
	  if (lastPos < 0) lastPos = fifoLength-1;
	  double valueSub = mapFifo[poli][lastPos]->GetBinContent(iPhi,iTheta); 
	  rollingMapAvg[rollingMapIndex(poli,iPhi,iTheta)] -= valueSub;
	  //	  if (poli==0 && iPhi==1 && iTheta==61 ) std::cout << "subtracting " << valueSub << std::endl;
	}	

	//and add the new addition to the fifo
	double valueAdd = mapFifo[poli][mapFifoPos]->GetBinContent(iPhi,iTheta);
	rollingMapAvg[rollingMapIndex(poli,iPhi,iTheta)] += valueAdd;
	//	if (poli==0 && iPhi==1 && iTheta==61 )  std::cout << "adding " << valueAdd << std::endl;

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
  
  //update flags and length
  noiseSummary->fifoLength = fifoLength;
  noiseSummary->mapFifoFillFlag = mapFifoFillFlag;
  noiseSummary->rmsFifoFillFlag = rmsFifoFillFlag;

  //update RMS stuff (always do this, it compresses well and is important)
  for (int phi=0; phi<NUM_PHI; phi++) {
    for (int ringi=0; ringi<NUM_ANTENNA_RINGS; ringi++) {
      for (int poli=0; poli<NUM_POLS; poli++) {
	AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
	AnitaPol::AnitaPol_t   pol  = (AnitaPol::AnitaPol_t)poli;
	noiseSummary->avgRMSNoise[phi][ringi][poli] = getAvgRMSNoise(phi,ring,pol);
      }
    }
  }


  //fill a map, but only if one of the flags is on, fillMap is used if both are selected
  if (fillAvgMap || fillMap) {
    for (int poli=0; (AnitaPol::AnitaPol_t)poli != AnitaPol::kNotAPol; poli++) {
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)poli;
      if (noiseSummary->avgMapProf[poli] != NULL) {
	delete noiseSummary->avgMapProf[poli];
	noiseSummary->avgMapProf[poli] = NULL;
      }
      //      std::cout << "filling map: " << poli << " " << mapFifoPos << " " << mapFifo[poli][mapFifoPos] << std::endl;
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


  for (int poli=0; poli<NUM_POLS; poli++) {
    for (int dir=0; dir<AnitaEventSummary::maxDirectionsPerPol; dir++) {
      double peakPhi = eventSummary->peak[poli][dir].phi;
      double peakTheta = eventSummary->peak[poli][dir].theta;
      
      if (mapFifo[0][0] == NULL) {//can't do it if you haven't filled anything yet
	eventSummary->peak[poli][dir].mapHistoryVal = -999;
      }
      else {
	int binPhi = mapFifo[poli][0]->GetXaxis()->FindBin(peakPhi);
	int binTheta = mapFifo[poli][0]->GetYaxis()->FindBin(-peakTheta);
	double avgNoise = rollingMapAvg[rollingMapIndex(poli,binPhi,binTheta)];
	eventSummary->peak[poli][dir].mapHistoryVal = avgNoise;
      }
    }
  }

  return;
}
