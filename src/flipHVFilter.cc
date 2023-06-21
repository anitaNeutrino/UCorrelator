#include "FFTtools.h" 
#include "FilterStrategy.h"
#include "flipHVFilter.h" 
#include "BasicFilters.h" 
#include "ResponseManager.h" 
#include "SystemResponse.h" 
#include <map> 
#include "TGraph.h" 
#include "AnalysisWaveform.h" 
#include "RawAnitaHeader.h"
#include "TRandom.h"
#include "TFile.h"
#include "Adu5Pat.h"
#include "DigitalFilter.h"
#include "RawAnitaHeader.h"
#include <math.h>// for isnan
#include "AntennaPositions.h"
#include <ctype.h>
#include <stdio.h>

UCorrelator::flipHVFilter::flipHVFilter(const char * index_file, bool anti)
  : anti(anti) 
{
  if (!index_file) 
  {
    TFile f(Form("%s/share/AnitaAnalysisFramework/responses/BH13TransferFn.root", getenv("ANITA_UTIL_INSTALL_DIR")));
    gPhase.push_back((TGraph*) f.Get("fixPhase"));
    gMag.push_back((TGraph*) f.Get("fixAmp"));
    f.Close();
  }

  ///load all the responses and the times for each time period
  else
  {
    TFile f(Form("%s/share/AnitaAnalysisFramework/responses/BH13TransferFn2020.root", getenv("ANITA_UTIL_INSTALL_DIR")));

    std::map<std::string,unsigned> notch_map; 


    TString index_fn; 
    index_fn.Form("%s/share/AnitaAnalysisFramework/responses/%s", getenv("ANITA_UTIL_INSTALL_DIR"),index_file); 

    std::ifstream inf(index_fn.Data());
    std::string notch_string; 
    unsigned tempTime; 
    if(inf)
    {
      while(inf >> notch_string >> tempTime)
      {
        end_times.push_back(tempTime); 

        if (!notch_map.count(notch_string))
        {
          //notch we haven't seen before 
          indices.push_back(gMag.size()); 
          notch_map[notch_string] = gMag.size(); 


          std::string mag_string = notch_string + "/mag"; 
          TGraph * mag = (TGraph*) f.Get(mag_string.c_str()); 
          gMag.push_back(mag); 

          std::string phase_string = notch_string + "/phase"; 
          TGraph * phase = (TGraph*) f.Get(phase_string.c_str()); 
          gPhase.push_back(phase); 
        }
        else 
        {
          indices.push_back(notch_map[notch_string]); 
        }
      }
    }
    else
    {
      fprintf(stderr,"Could not open file %s\n", index_fn.Data());
    }
  }
}


UCorrelator::flipHVFilter::~flipHVFilter() {

  for (unsigned i = 0; i < gPhase.size(); i++) {
  
    delete gPhase[i];
    delete gMag[i];
  }
}


void UCorrelator::flipHVFilter::process(FilteredAnitaEvent * ev) {

  processOne(getWf(ev, 44, AnitaPol::kVertical), ev->getHeader(), 44, AnitaPol::kVertical);

}


void UCorrelator::flipHVFilter::processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol) {

  if(whichAnt != 44 || whichPol != 0) return;

  int index = 0; 

  //check curent index 
  if (end_times.size()) index = indices[std::upper_bound(end_times.begin(), end_times.end(), header->triggerTime) - end_times.begin()];

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  int old_size = awf->Neven();
  int nf = awf->Nfreq();
  double df = awf->deltaF();
  
  for( int i =0; i < nf; i++) {
	
    double f =i*df; 
    double lookup_f = f;
    
    if (!end_times.size()) lookup_f*=1e9; // old version used Hz, not GHz
    
    double phase = anti ? awf->freq()[i].getPhase()-FFTtools::evalEvenGraph(gPhase[index], lookup_f) : awf->freq()[i].getPhase()+FFTtools::evalEvenGraph(gPhase[index],lookup_f);

    double mag = anti ? awf->freq()[i].getAbs()*FFTtools::evalEvenGraph(gMag[index],lookup_f) : awf->freq()[i].getAbs()*(1./gMag[index]->Eval(lookup_f));
    
    if( f>=.1 && f<=1.3) awf->updateFreq()[i].setMagPhase(mag, phase);

    }
    
    awf->updateEven()->Set(old_size);
}
