#ifndef _UCORRELATOR_BASELINE_H 
#define _UCORRELATOR_BASELINE_H 

/* Implementation of RF Baselines */ 

class TGraph; 
#include "AnitaConventions.h" 
#include "AnitaGeomTool.h"
namespace UCorrelator
{


  class Baseline
  {

  /** This class is used to build, save and fetch baselines
   * for use with the adaptive filter. 
   *  
   *  Baselines are per run and appear to be simply averaged averforms. 
   *
   * */ 


    public: 
      /** Create a baseline for the given run using the given number of traces
       * to average
       *
       * If persistdir is given, that directory will first be checked for an
       * already computed baseline with navg and that will be loaded
       * instead, if it exists.  If a baseline is not found in persistdir, the
       * baseline will be generated and saved there. 
       *
       * The environmental variable ANITA_ROOT_DATA must be defined for this to work properly. 
       *
       *
       * The baselines are averaged magnitudes (not powers) in "raw units", not dB (and not normalized by the number of samples)
       */ 

      Baseline(int run,  int navg = 5000, 
                const char * persistdir = ""); 

      /** Save baseline to directory. If persistdir is given in the constructor,
       * the baseline will already have been saved there so you probably 
       * don't need to call this .
       *
       * */ 
      void saveToDir(const char * dir); 

      /** Getters for the baseline */ 
      const TGraph *getBaseline(AnitaPol::AnitaPol_t pol, int ant)
        { return pol == AnitaPol::kHorizontal? hpol[ant] : vpol[ant] ; }
      const TGraph *getBaselineHPol(int ant)  { return hpol[ant]; } 
      const TGraph *getBaselineVPol(int ant)  { return vpol[ant]; } 

      const TGraph *getBaselineAverage(AnitaPol::AnitaPol_t pol)
        { return pol == AnitaPol::kHorizontal? hpol_avg : vpol_avg ; }
      const TGraph *getBaselineAverageHPol()  { return hpol_avg; } 
      const TGraph *getBaselineAverageVPol()  { return vpol_avg; } 

    private: 
      TGraph * hpol[NUM_SEAVEYS]; 
      TGraph * vpol[NUM_SEAVEYS];  
      TGraph * hpol_avg;
      TGraph * vpol_avg;
      int navg; 
      int run; 


  };

}

#endif 
