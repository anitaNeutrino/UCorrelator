#ifndef UCORRELATOR_CORRELATOR_H
#define UCORRELATOR_CORRELATOR_H


#include "AnitaConventions.h"
#include "TH2.h"

class FilteredAnitaEvent; 
class AnalysisWaveform; 

#define NANTENNAS (NUM_DIGITZED_CHANNELS/2)

namespace UCorrelator
{

  class Correlator
  {
    public:
      Correlator(int nphi, double phimin, double phimax, int ntheta, double theta_min, double theta_max); 
      void compute(const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t pol); 
      const TH2 * getHist() const { return &hist; }  
      
      TH2D* computeZoomed(double phi, double theta, int ntheta, double dtheta, int nphi, double dphi, int nant, TH2D * useme = 0); 

      void setDisallowedAntennas(uint64_t disallowed) { disallowed_antennas = disallowed; } 
      void setAllowedAntennas(uint64_t allowed) { disallowed_antennas = ~allowed; } 
      void setMaxAntennaMaxPhiDistance(double max_ant_phi) { max_phi = max_ant_phi; } 

      void setGroupDelayFlag(bool flag) { groupDelayFlag = flag; } 

      void setPadFactor(int pad) { pad_factor = pad; } 
      virtual ~Correlator(); 

    private: 
      AnalysisWaveform* padded_waveforms[NANTENNAS]; 
      AnalysisWaveform* correlations[NANTENNAS][NANTENNAS]; 
      double rms[NANTENNAS]; 
      TH2D hist; 
      TH2I norm; 

      double max_phi;
      uint64_t disallowed_antennas;
      double * phi_table; 
      double * theta_table; 
      double * cos_phi_table; 
      double * tan_theta_table; 
      double * cos_theta_table; 
      int pad_factor;
      const FilteredAnitaEvent * ev; 
      AnitaPol::AnitaPol_t pol; 
      bool groupDelayFlag; 

      AnalysisWaveform * getCorrelation(int ant1, int ant2); 
      void compute_table(); 
      void doAntennas(int ant1, int ant2, TH2D * hist, TH2I * norm, bool useArray); 
      void reset(); 
      double getDeltaTFast(int ant1, int ant2, int phibin, int thetabin); 
      double getDeltaT(int ant1, int ant2, double phi, double theta); 
      double getGroupDelay(double phidiff, double theta); 
  }; 



}



#endif
