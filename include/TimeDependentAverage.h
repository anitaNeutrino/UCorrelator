#ifndef _UCORRELATOR_TIME_DEPENDENT_AVERAGE_H 
#define _UCORRELATOR_TIME_DEPENDENT_AVERAGE_H 

/* Implementation of RF TimeDependentAverages */ 

class TH2F; 
class TH1D; 
class TH1I; 
class TH1; 
class TH2D; 
#include "AnitaConventions.h" 
#include "AnitaGeomTool.h"
#include "TMutex.h" 
  
namespace UCorrelator
{

  /** This class is used to build, save and fetch spectrum averages
   * for use with the adaptive filter and other time dependent stuff. 
   *  
   *
   *  It also calculates something I call the peakiness.  
   *
   * */ 



  class TimeDependentAverage
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


      TimeDependentAverage(int run, int nsecs = 10, 
                      const char * persistdir = 0, 
                      double max_bottom_top_ratio =4.0, int min_norm = 5, double max_power = 1e6); 


      static const TimeDependentAverage * defaultThermal(); 
      virtual ~TimeDependentAverage(); 

      /** compute the peakiness with the given thermal spectrum average and fraction of spectrum used for normalization
       * If 0 is passed for thermal, a default thermal spectrum  is used from calibration terminated amp data. 
       * */ 
      void computePeakiness(const TimeDependentAverage * thermal = 0 , double fractionForNormalization = 0.5) const;  

      /** Save average to directory. If persistdir is given in the constructor,
       * the baseline will already have been saved there so you probably 
       * don't need to call this .
       *
       * */ 
      void saveToDir(const char * dir); 

      /* Get the spectrogrma at a time. This currently does not interpoolate... just picks the closet bin */ 
      TH1 *getSpectrumAverage(AnitaPol::AnitaPol_t pol, int ant, double t, bool db = false, bool minbias = false) const;

      /* Get the spectrum at a given percentile. e.g. at pct = 0.5, gives the median spectrum */ 
      TH1 *getSpectrumPercentile(AnitaPol::AnitaPol_t pol, int ant, double pct = 0.1, bool db = false, bool minbias = false) const; 

      double getStartTime() const; 
      double getEndTime() const; 

      const TH2F * getSpectrogram(AnitaPol::AnitaPol_t pol, int ant, bool minbias = false) const { return minbias ? avgs_minbias[ant][pol] : avgs[ant][pol] ; }
      const TH1D * getRMS(AnitaPol::AnitaPol_t pol, int ant) const { return rms[ant][pol] ; } 
      double getRMS(AnitaPol::AnitaPol_t pol, int ant, double t) const; 
      const TH1I * getNBlasts() const { return nblasts; } 
      double getBlastFraction(double t) const; 
      const TH1I * getNorms(bool minbias = false) const { return minbias ? norms_minbias : norms; }
      const TH2D * getPeakiness(AnitaPol::AnitaPol_t pol, int ant, bool minbias = false) const; 
      int getRun() const { return run; } 
      int getNsecs() const { return nsecs; } 

    private: 
      TH2F * avgs[NUM_SEAVEYS][2]; 
      TH2F * avgs_minbias[NUM_SEAVEYS][2]; 
      TH1D * rms[NUM_SEAVEYS][2]; 
      TH1I * nblasts; 
      TH1I * norms; 
      TH1I * norms_minbias; 

      mutable TMutex m; 
      mutable TH2D * peakiness[NUM_SEAVEYS][2]; 
      mutable TH2D * peakiness_minbias[NUM_SEAVEYS][2]; 
      int computeAverage(double max_r, int min_norm, double max_power); 
      int nsecs; 
      int run; 
      const char * dir; 

  };


  /* This will try to to find an appropriate spectrum average
   * for a given time. Useful for MC or other cases where we span runs */ 

  class TimeDependentAverageLoader
  {

    public: 
      TimeDependentAverageLoader(const char * dir = 0, int nsecs = 10); 

      /* Tries to find a TimeDependentAverage with the right run */ 
      const TimeDependentAverage * avg(double t) const; 


      /** Static member functions, use a time dependent average loader with the environmental variable  in the background if NULL passed
       *
       * These are thread safe, but if you change nsecs in between cals, it will be very inefficient. 
       * */ 
      static double getRMS(double t, AnitaPol::AnitaPol_t pol, int ant, int nsecs = 10);  
      static double getPayloadBlastFraction(double t,  int nsecs =10);  
      
      int getNsecs() const { return nsecs; } 
    private: 
      mutable TimeDependentAverage * tavg; 
      const char * dir; 
      int nsecs; 
  }; 


}

#endif 
