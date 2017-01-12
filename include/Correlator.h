#ifndef UCORRELATOR_CORRELATOR_H
#define UCORRELATOR_CORRELATOR_H


#include "AnitaConventions.h"
#include "TH2.h"
#include <stdint.h>

class FilteredAnitaEvent; 
class AnalysisWaveform; 
class TrigCache; 

#define NANTENNAS NUM_SEAVEYS

namespace UCorrelator
{
  
  //This is just to hide the OpenMP stuff from cling
  class CorrelatorLocks; 

  /** This creates the inteferometric map for an ANITA event */ 
  class Correlator
  {
    public:
      /** Create a correlator with the following options for the rough map */ 
      Correlator(int nphi, double phimin, double phimax, int ntheta, double theta_lowest, double theta_highest, bool use_bin_center = false, bool scale_by_cos_theta = false, double baseline_weight = 0); 

      /** Compute the rough correlation map for the event and pol */ 
      void compute(const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t pol); 

      /** Get the rough correlation map */ 
      const TH2D * getHist() const { return hist; } 
 
      /** Get the rough correlation map normalization */ 
      const TH2I * getNorm() const { return norm; } 
      

      /** Compute a zoomed in map around phi and theta. nphi,dphi,ntheta,dtheta. If nant is non-zero, only the nearest nant antennas are used. You can use useme to avoid allocating a new TH2.  */ 
      TH2D* computeZoomed(double phi, double theta, int nphi, double dphi,  int ntheta, double dtheta, int nant = 0, TH2D * useme = 0); 

      /** Disable the antennas given by the bitmap */ 
      void setDisallowedAntennas(uint64_t disallowed) { disallowed_antennas = disallowed; } 

      /** Enable only the antennas given by the bitmap */ 
      void setAllowedAntennas(uint64_t allowed) { disallowed_antennas = ~allowed; } 

      /** An antenna only contributes to an angle if it's within max_phi  of it */
      void setMaxAntennaMaxPhiDistance(double max_ant_phi) { max_phi = max_ant_phi;  max_phi2 = max_phi * max_phi; } 

      /** Use the group delay in computing the delay */ 
      void setGroupDelayFlag(bool flag) { groupDelayFlag = flag; } 


      /** Get the correlation between two antennas */
      const AnalysisWaveform * getCorrelationGraph(int ant1, int ant2) { return getCorrelation(ant1,ant2); }

      /** Set the supersampling factor */ 
      void setPadFactor(int pad) { pad_factor = pad; } 

      /** Debugging method to dump out some info to a file */ 
      void dumpDeltaTs(const char * file) const; 
      virtual ~Correlator(); 

    private: 
      AnalysisWaveform* padded_waveforms[NANTENNAS]; 
      AnalysisWaveform* correlations[NANTENNAS][NANTENNAS]; 

      TH2D *hist; //wow, apparently sizeof(TH2D) is huge... that's why this is on the heap 
      TH2I *norm; 

      TrigCache * trigcache[NUM_ANITAS+1]; 
      double rms[NANTENNAS]; 

      double max_phi, max_phi2;
      uint64_t disallowed_antennas;
      int pad_factor;
      const FilteredAnitaEvent * ev; 
      AnitaPol::AnitaPol_t pol; 
      bool groupDelayFlag; 
      bool use_bin_center; 
      bool scale_cos_theta; 
      double baselineWeight;

      AnalysisWaveform * getCorrelation(int ant1, int ant2); 
      void doAntennas(int ant1, int ant2, TH2D * hist, TH2I * norm, const TrigCache * tc, const double * center_point  = 0); 
      void reset(); 

      CorrelatorLocks * locks; 
  }; 



}



#endif
