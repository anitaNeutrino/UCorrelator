#ifndef _UCORRELATOR_SPECTRUM_AVERAGE_H 
#define _UCORRELATOR_SPECTRUM_AVERAGE_H 

/* Implementation of RF SpectrumAverages */ 

class TH2D; 
class TH1; 
#include "AnitaConventions.h" 
#include "AnitaGeomTool.h"
  
namespace UCorrelator
{

  /** This class is used to build, save and fetch spectrum averages
   * for use with the adaptive filter. 
   *  
   *
   *  It also calculates something I call the peakiness.  
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
                      const char * persistdir = 0,const char * selection=0,
                      double max_bottom_top_ratio =1.5); 


      static const SpectrumAverage * defaultThermal(); 
      virtual ~SpectrumAverage(); 

      /** compute the peakiness with the given thermal spectrum average and fraction of spectrum used for normalization
       * If 0 is passed for thermal, a default thermal spectrum  is used from calibration terminated amp data. 
       * */ 
      void computePeakiness(const SpectrumAverage * thermal = 0 , double fractionForNormalization = 0.5);  

      /** Save average to directory. If persistdir is given in the constructor,
       * the baseline will already have been saved there so you probably 
       * don't need to call this .
       *
       * */ 
      void saveToDir(const char * dir); 

      TH1 *getSpectrumAverage(AnitaPol::AnitaPol_t pol, int ant, double t, bool db = false) const;
      TH1 *getSpectrumPercentile(AnitaPol::AnitaPol_t pol, int ant, double pct = 0.1, bool db = false) const; 

      double getStartTime() const; 
      double getEndTime() const; 

      const TH2D * getSpectrogram(AnitaPol::AnitaPol_t pol, int ant) const {return avgs[ant][pol]; } 
      const TH2D * getPeakiness(AnitaPol::AnitaPol_t pol, int ant) const { return peakiness[ant][pol]; } 
      int getRun() const { return run; } 
      int getNsecs() const { return nsecs; } 

    private: 
      TH2D * avgs[NUM_SEAVEYS][2]; 
      TH2D * peakiness[NUM_SEAVEYS][2]; 
      int computeAverage(const char * selection, double max_r); 
      int nsecs; 
      int run; 

  };

}

#endif 
