#include <iostream>
#include "FFTtools.h"
#include "AnitaEventSummary.h"
#include "UCUtil.h"
#include "TTree.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "RawAnitaHeader.h"

using namespace std;

void doSummaryResolutions(int firstRun, int lastRun, const char * folder)
{

 const char * filter = "sinsub_5_3_ad_2";
  
  AnitaPol::AnitaPol_t pol;

  char cpol[100];
  int whichPulser = 0;  
  std::string pulser="";

  if (firstRun>300){ // WAIS
    pol = AnitaPol::kHorizontal;
    pulser+="WAIS";
  } else if (firstRun<150){ // LDB VPOL 
    pol = AnitaPol::kVertical;
    whichPulser = 1;
    pulser+="LDB";
  } else if (firstRun<154){ // LDB HPOL
    pol = AnitaPol::kHorizontal;
    whichPulser = 1;
    pulser+="LDB";    
  } else if (firstRun<172){ // LDB VPOL
    pol = AnitaPol::kVertical;
    whichPulser = 1;
    pulser+="LDB";
  } else if (firstRun>203 && firstRun<251){ // QUIET TIME TO STUDY THE SUN
    pol = AnitaPol::kHorizontal;
    whichPulser = 2;
    pulser+="SUN";
  } else {
    std::cout << "Unknown firstRun" << std::endl;
    return;
  }

  if (pol == AnitaPol::kVertical) sprintf(cpol, "VPOL");
  else if (pol == AnitaPol::kHorizontal) sprintf(cpol, "HPOL");

  TChain *chain = new TChain("resolution");
  
  for (int irun=firstRun; irun<=lastRun;irun++){
    TString inname; 
    inname.Form("%s/Resolution_%s_%s_%d_%s.root", folder, pulser.c_str(), cpol, irun, filter);
    cout << inname << endl;
    chain->Add(inname);
  }

  AnitaEventSummary * sum = new AnitaEventSummary;   
  chain->SetBranchAddress("summary",    &sum);

  Long_t nEntries=chain->GetEntries();
  cout << "There are " << nEntries << " in this chain" << endl;
  chain->GetEntry(0);
  Double_t firstTS = sum->eventNumber;//(double)triggerTime;
  chain->GetEntry(nEntries-1);
  Double_t lastTS = sum->eventNumber;//(double)triggerTime;
  cout << "From time " << firstTS << " " << " to time " << lastTS << endl;
  
  string angleNames[2] = {"phi", "theta"};
  double angleMin[2]   = {  -5.,     -2.};
  double angleMax[2]   = {   5.,      2.};
  
  string otherVarsName[4] = {"phi", "theta", "snr", "eventNumber"};
  double otherVarMin[4]   = {   0.,      0.,    0.,    firstTS};
  double otherVarMax[4]   = { 360.,     10.,   20.,     lastTS};
  
  // array of 1D histograms with the delta phi or delta theta
  TH1D *hDelta[2];

  // array of 2D histograms:
  // y axis is the delta phi/theta
  // x axis is the other variable
  TH2D *hDelta2D[2][4];
  
  for (int iangle=0; iangle<2; iangle++){
    hDelta[iangle] = new TH1D(Form("hDelta_%i", iangle), "",  100, angleMin[iangle], angleMax[iangle]);

    for (int ivar=0; ivar<4; ivar++){
      hDelta2D[iangle][ivar] = new TH2D(Form("hDelta2D_%i_%i", iangle, ivar), "", 100, otherVarMin[ivar], otherVarMax[ivar], 100, angleMin[iangle], angleMax[iangle]);

    }
    
  }


  double anglesMeas[2];
  double anglesTrue[2];
  double otherVars[4];
  double delta;
  
  for (int ientry=0; ientry<nEntries; ientry++){
    chain->GetEntry(ientry);

    anglesMeas[0] = sum->peak[pol][0].phi;
    anglesMeas[1] = sum->peak[pol][0].theta;

    if (whichPulser==2){
      anglesTrue[0] = sum->sun.phi;
      anglesTrue[1] = sum->sun.theta;
    } if (whichPulser==1){
      anglesTrue[0] = sum->ldb.phi;
      anglesTrue[1] = sum->ldb.theta;
    }else{
      anglesTrue[0] = sum->wais.phi;
      anglesTrue[1] = sum->wais.theta;      
    }

    otherVars[0] = anglesMeas[0];
    otherVars[1] = anglesMeas[1];
    otherVars[2] = sum->coherent[pol][0].snr;
    otherVars[3] = (double)sum->eventNumber;
    
    
    for (int iangle=0; iangle<2; iangle++){
      delta = anglesMeas[iangle] - anglesTrue[iangle];
      //      cout << anglesMeas[iangle] << " " <<  anglesTrue[iangle] << " " << delta << endl;
      hDelta[iangle]->Fill(delta);
      for (int ivar=0; ivar<4; ivar++){
	hDelta2D[iangle][ivar]->Fill(otherVars[ivar], delta);
      }
    }
    

  }

  TFile *fout = new TFile(Form("%s/SummaryResolutions_%s.root", folder, pulser.c_str()), "recreate");
  
    
  for (int iangle=0; iangle<2; iangle++){
    hDelta[iangle]->Write(("hDelta_"+angleNames[iangle]).c_str());
    for (int ivar=0; ivar<4; ivar++){
      hDelta2D[iangle][ivar]->Write(("hDelta2d_"+angleNames[iangle]+"_"+otherVarsName[ivar]).c_str());
    }
  }
  
}


int main (int nargs, char ** args)
{
   
  int firstRun        = nargs < 2 ? 345 : atoi(args[1]); 
  int lastRun         = nargs < 3 ? 345 : atoi(args[2]); 
  const char * folder = nargs < 4 ? 0 : args[3]; 
  
  doSummaryResolutions(firstRun, lastRun, folder); 


}
