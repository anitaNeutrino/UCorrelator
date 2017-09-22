#include "SineSubtractCache.h"
#include <iostream>
#include "AnitaDataset.h"
#include "TFile.h"
#include "TTree.h"
#include "FilteredAnitaEvent.h"
#include "SineSubtractCache.h"
#include "FilterStrategy.h"
#include "RawAnitaHeader.h"
#include "TDirectory.h"


const TString ssrTreeName = "sineSubResultTree";

TString UCorrelator::SineSubtractCache::branchName(AnitaPol::AnitaPol_t pol, Int_t ant){
  return TString::Format("ssr_%d_%d", (int)pol, ant);
}

TString UCorrelator::SineSubtractCache::fileName(const char* specDir, UInt_t hash, Int_t run){
  int av = AnitaVersion::get();
  return TString::Format("%s/sineSubResults_%u_anita%d_run%d.root", specDir, hash, av, run);
}


void UCorrelator::SineSubtractCache::makeCache(int run, SineSubtractFilter* ssf){

  UCorrelator::SineSubtractFilter::setUseCache(false);
  FFTtools::SineSubtractResult* results[AnitaPol::kNotAPol][NUM_SEAVEYS] = {{NULL}};
  
  const char* ssDesc = ssf->description();
  if(!ssDesc){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't cache or use cached results without SineSubtract description!" << std::endl;
  }  
  else{
    UInt_t hash = TString(ssDesc).Hash();
    
    const char* specDir = getenv("UCORRELATOR_SPECAVG_DIR");
    if(!specDir){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't cache or use cached results without UCORRELATOR_SPECAVG_DIR environment variable!" << std::endl;
    }
    else{

      std::cout << "Info in " << __PRETTY_FUNCTION__ << ", will generate cached sine subtraction results!" << std::endl;
      AnitaDataset d(run);

      FilterStrategy fs;
      fs.addOperation(ssf);
  
      TFile* fOut = new TFile(fileName(specDir, hash, run), "recreate");
      TTree* tOut = new TTree(ssrTreeName, ssrTreeName);
      UInt_t eventNumber = 0;
      tOut->Branch("eventNumber", &eventNumber);
      for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
        AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
        for(int ant=0; ant < NUM_SEAVEYS; ant++){
          tOut->Branch(branchName(pol, ant), &results[pol][ant]);
        }
      }


      const int n = d.N();
      const double deltaPrint = double(n)/1000;
      double nextPrint = 0;
      for(int entry=0; entry < d.N(); entry++){
        d.getEntry(entry);

        FilteredAnitaEvent fEv(d.useful(), &fs, d.gps(), d.header());

        for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
          AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;          
          for(int ant=0; ant < NUM_SEAVEYS; ant++){
            const FFTtools::SineSubtract* ss = ssf->sinsub(pol, ant);
            FFTtools::SineSubtractResult* ssr = const_cast<FFTtools::SineSubtractResult*>(ss->getResult());
            results[pol][ant] = ssr;
          }
        }

        eventNumber = d.header()->eventNumber;
        tOut->Fill();
        if(entry >= nextPrint){
          const int nm = 50;
          int m = nm*nextPrint/n;
          fprintf(stderr, "\r%4.2f %% complete", 100*nextPrint/n);
          std::cerr << "[";
          for(int i=0; i < m; i++) std::cerr << "=";
          for(int i=m; i < nm; i++) std::cerr << " ";
          std::cerr << "]";
          nextPrint += deltaPrint;
        }
        // if(entry > 100) break;
      }
      tOut->BuildIndex("eventNumber"); // does this get saved?
      fOut->Write();
      fOut->Close();

    }
  }
}









/** 
 * Constructor
 * 
 * @param ssDesc description string, used for identifying file of cached results
 */
UCorrelator::SineSubtractCache::SineSubtractCache(const char* ssDesc)
    : fFile(NULL), fTree(NULL), fDescHash(0), fSpecDir(), fCurrentRun(-1), fLastEventNumber(0){


  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      results[pol][ant] = NULL;
    }
  }
  
  if(!ssDesc){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't cache or use cached results without SineSubtract description!" << std::endl;
  }  
  else{
    fDescHash = TString(ssDesc).Hash();
    
    const char* specDir = getenv("UCORRELATOR_SPECAVG_DIR");
    if(!specDir){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't cache or use cached results without UCORRELATOR_SPECAVG_DIR environment variable!" << std::endl;
    }
    else{
      fSpecDir = specDir;
    }
  }
}



UCorrelator::SineSubtractCache::~SineSubtractCache(){

  if(fFile){
    fFile->Close();
  }
}



const FFTtools::SineSubtractResult* UCorrelator::SineSubtractCache::getResult(UInt_t eventNumber, AnitaPol::AnitaPol_t pol, Int_t ant){

  // hard to check whether anita version is correct...
  // this should happen

  if(eventNumber != fLastEventNumber){
    int run = AnitaDataset::getRunContainingEventNumber(eventNumber);
    if(run!=fCurrentRun){
      loadRun(run);
    }
    if(fTree){
      Int_t entry = fTree->GetEntryNumberWithIndex(eventNumber);
      std::cerr << entry << std::endl;
      if(entry >= 0){
        fTree->GetEntry(entry);

        std::cerr << fCurrentRun << "\t" << fLastEventNumber << std::endl;
        
      }
      else{
        std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't find entry "
                  << entry << " in " << fTree->GetName() << " in file " << fFile->GetName()
                  << " for eventNumber " << eventNumber << std::endl;
        return NULL;
      }
    }
  }
  
  return results[pol][ant];
}



void UCorrelator::SineSubtractCache::loadRun(Int_t run){

  if(fSpecDir.Length() > 0 && fDescHash > 0){
    if(fFile){
      fFile->Close();
      fFile = NULL;
      fTree = NULL;
    }

    const TString theRootPwd = gDirectory->GetPath();
    TString fName = fileName(fSpecDir, fDescHash, run);
    fFile = TFile::Open(fName, "read");
    if(fFile){
      fTree = (TTree*) fFile->Get(ssrTreeName);

      fTree->SetBranchAddress("eventNumber", &fLastEventNumber);
      for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
        AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
        for(int ant=0; ant < NUM_SEAVEYS; ant++){
          fTree->SetBranchAddress(branchName(pol, ant), &results[pol][ant]);
        }
      }
      fTree->BuildIndex("eventNumber");
      fTree->GetEntry(0);
      fCurrentRun = run;
      
      std::cerr << fCurrentRun << "\t" << fLastEventNumber << std::endl;
      

      gDirectory->cd(theRootPwd); 
    }
    else{
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't open file " << fName << std::endl;
    }
  }
}
