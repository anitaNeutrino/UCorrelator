#ifndef UCORRELATOR_CORRELATOR_H
#define UCORRELATOR_CORRELATOR_H


#include "AnitaConventions.h"
#include "TH2.h"

class FilteredAnitaEvent; 
class AnalysisWaveform; 
class TrigCache; 

#define NANTENNAS NUM_SEAVEYS

namespace UCorrelator
{

  class Correlator
  {
    public:
      Correlator(int nphi, double phimin, double phimax, int ntheta, double theta_min, double theta_max); 
      void compute(const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t pol); 
      const TH2 * getHist() const { return &hist; }  
      const TH2 * getNorm() const { return &norm; } 
      
      TH2D* computeZoomed(double phi, double theta, int nphi, double dphi,  int ntheta, double dtheta, int nant, TH2D * useme = 0); 

      void setDisallowedAntennas(uint64_t disallowed) { disallowed_antennas = disallowed; } 
      void setAllowedAntennas(uint64_t allowed) { disallowed_antennas = ~allowed; } 
      void setMaxAntennaMaxPhiDistance(double max_ant_phi) { max_phi = max_ant_phi; } 

      void setGroupDelayFlag(bool flag) { groupDelayFlag = flag; } 
      const AnalysisWaveform * getCorrelationGraph(int ant1, int ant2) { return getCorrelation(ant1,ant2); }

      void setPadFactor(int pad) { pad_factor = pad; } 
      void dumpDeltaTs(const char * file) const; 
      virtual ~Correlator(); 

    private: 
      AnalysisWaveform* padded_waveforms[NANTENNAS]; 
      AnalysisWaveform* correlations[NANTENNAS][NANTENNAS]; 
      TrigCache * trigcache; 
      double rms[NANTENNAS]; 
      TH2D hist; 
      TH2I norm; 

      double max_phi;
      uint64_t disallowed_antennas;
      int pad_factor;
      const FilteredAnitaEvent * ev; 
      AnitaPol::AnitaPol_t pol; 
      bool groupDelayFlag; 

      AnalysisWaveform * getCorrelation(int ant1, int ant2); 
      void doAntennas(int ant1, int ant2, TH2D * hist, TH2I * norm, const TrigCache * tc); 
      void reset(); 

  }; 



}



#endif
