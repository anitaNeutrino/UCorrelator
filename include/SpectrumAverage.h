#ifndef _UCORRELATOR_SPECTRUM_AVERAGE_H 
#define _UCORRELATOR_SPECTRUM_AVERAGE_H 

/* Implementation of RF SpectrumAvereages */ 

class TH2; 
#include "AnitaConventions.h" 
#include "AnitaGeomTool.h"
namespace UCorrelator
{

  /** This class is used to build, save and fetch spectrum averages
   * for use with the adaptive filter. 
   *  
   * */ 


  class SpectrumAverage
  {


    public: 
      /** Create a time-dependent average
       *
       * If persistdir is given, that directory will first be checked for an
       * If a baseline is not found in persistdir, the
       * baseline will be generated and saved there. 
       *
       * The environmental variable ANITA_ROOT_DATA must be defined for this to work properly. 
       *
       */ 


      SpectrumAverage(int run, int nsecs = 60, 
                      const char * persistdir = "",const char * selection="",
                      double max_bottom_top_ratio =1.5); 


      /** Save average to directory. If persistdir is given in the constructor,
       * the baseline will already have been saved there so you probably 
       * don't need to call this .
       *
       * */ 
      void saveToDir(const char * dir); 

      TH2 *getSpectrumAverage(AnitaPol::AnitaPol_t pol, int ant, double t);  
      TH2 *getSpectrumPercentile(AnitaPol::AnitaPol_t pol, int ant, double pct = 0.1); 
      const TH2 * getSpectogram() const; 

    private: 
      TH2 * avgs[NUM_SEAVEYS][2]; 
      int nsecs; 
      int run; 

  };

}

#endif 
