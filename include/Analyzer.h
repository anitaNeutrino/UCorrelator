#ifndef UCORRELATOR_ANALYZER_H
#define UCORRELATOR_ANALYZER_H

#include "Correlator.h" 
#include "WaveformCombiner.h"
#include <vector>
#include "TH2D.h"
#include "simpleStructs.h"

#include "AnitaEventSummary.h" 
#include "ResponseManager.h" 

class FilteredAnitaEvent; 
class UsefulAdu5Pat; 
class AnalysisWaveform; 
class TPad; 

namespace FFTtools
{
  class DigitalFilter; 
}

namespace UCorrelator
{

  class AnalysisConfig; 

  /** 
   * The analyzer is the main workhorse of UCorrelator. It plumbs together various other parts of UCorrelator to do stuff.  
   *
   */ 

  class Analyzer
  {

    public:
      /** Create an Analyzer. By default, the default AnalysisConfig parameters are used, but one may be passed. If interactive is passed, 
       * more internal state will be accessible */ 
      Analyzer(const AnalysisConfig * cfg = 0, bool interactive = false); 

      /** Destructor */ 
      virtual ~Analyzer();



      /** Analyze the event and put the results into the summary. The summary is reinitialized, clobbering everything. */
      void analyze(const FilteredAnitaEvent * event, AnitaEventSummary *summary); 

      /** Retrieve the internal correlator. Note that if multiple polarizations are analyzed, the Correlator's internal state
       *  will be related to the last polarization used. */ 
      Correlator * getCorrelator() { return &corr; } 

      /** Retrieve the internal waveform combiner. Note that if multiple polarizations are analyzed, the WaveformCombiner's internal state
       *  will be related to the last polarization used. */ 
      WaveformCombiner * getWaveformCombiner() { return &wfcomb; } 

      /** Retrieve the internal xpol waveform combiner. Note that if multiple polarizations are analyzed, the WaveformCombiner's internal state
       *  will be related to the last polarization used. */ 
      WaveformCombiner * getXPolWaveformCombiner() { return &wfcomb_xpol; } 


      /** Return the correlation map for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const TH2 * getCorrelationMap(AnitaPol::AnitaPol_t pol) const  { return correlation_maps[pol] ; } 

      /** Return the ith zoomed correlation map for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const TH2 * getZoomedCorrelationMap(AnitaPol::AnitaPol_t pol, int i) const { return zoomed_correlation_maps[pol][i]; }

      /** Return the ith coherent waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getCoherent(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return coherent[pol][i][filtered ? 1 : 0]; } 

       /** Return the ith deconvolved waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getDeconvolved(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return deconvolved[pol][i][filtered ? 1: 0]; } 

      /** Return the ith coherent xpol waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getCoherentXpol(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return coherent_xpol[pol][i][filtered ? 1: 0]; } 

       /** Return the ith deconvolved xpol waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getDeconvolvedXpol(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return deconvolved_xpol[pol][i][filtered ? 1: 0]; } 

      /** Return the ith coherent averaged power for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const TGraphAligned * getCoherentPower(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return coherent_power[pol][i][filtered ? 1: 0]; } 

       /** Return the ith deconvolved averaged power for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const TGraphAligned * getDeconvolvedPower(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return deconvolved_power[pol][i][filtered ? 1: 0]; } 

      ResponseManager * getResponseManager() { return &responses; } 




      /** Populate the pads with a summary of the pointing. Only makes sense if interactive mode is on. If the analyzer wasn't instructed to do 
       * both polarities, then it won't populate any it wasn't instructed to do. If 0 or NULL is passed, a new canvas is made. */ 
      void drawSummary(TPad *chpol = 0, TPad * cvpol = 0) const; 


      double getRoughPhi(AnitaPol::AnitaPol_t pol, int i) const { return rough_peaks[pol][i].first; }
      double getRoughTheta(AnitaPol::AnitaPol_t pol, int i) const { return -rough_peaks[pol][i].second; }

      void clearInteractiveMemory(double frac = 0.5) const; 
    private:

      void fillWaveformInfo(const AnalysisWaveform * wf, const AnalysisWaveform * xpol_wf, const TGraph * power, AnitaEventSummary::WaveformInfo * info, AnitaPol::AnitaPol_t pol); 
      void fillPointingInfo(double rough_phi, double rough_theta, AnitaEventSummary::PointingHypothesis * point,
                            UsefulAdu5Pat * pat, double hwAngle, UShort_t triggered_sectors, UShort_t masked_sectors, UShort_t triggered_sectors_xpol, UShort_t masked_sectors_xpol); 
      void fillFlags(const FilteredAnitaEvent * fae, AnitaEventSummary::EventFlags * flags, UsefulAdu5Pat * pat); 

      TH2D* correlation_maps[2]; 
      std::vector<TH2D*>  zoomed_correlation_maps[2]; 
      std::vector<AnalysisWaveform *> coherent[2][2]; 
      std::vector<AnalysisWaveform *> deconvolved[2][2]; 
      std::vector<TGraphAligned *> coherent_power[2][2]; 
      std::vector<TGraphAligned *> deconvolved_power[2][2]; 
      std::vector<AnalysisWaveform *> coherent_xpol[2][2]; 
      std::vector<AnalysisWaveform *> deconvolved_xpol[2][2]; 
      std::vector<TGraphAligned *> coherent_power_xpol[2][2]; 
      std::vector<TGraphAligned *> deconvolved_power_xpol[2][2]; 
      std::vector<std::pair<double,double> > rough_peaks[2]; 
      AnitaEventSummary last; //used in interactive mode by drawSummary
      TGraph*  avg_spectra[2]; 

      const AnalysisConfig * cfg; 
      Correlator corr; 
      ResponseManager responses; 
      WaveformCombiner wfcomb; 
      WaveformCombiner wfcomb_xpol; 
      WaveformCombiner wfcomb_filtered; 
      WaveformCombiner wfcomb_xpol_filtered; 
      TH2D *zoomed; 
      double maprms; 
      FFTtools::DigitalFilter *  power_filter; 
      bool interactive; 
      bool interactive_deconvolved; 
      bool interactive_xpol_deconvolved; 

      mutable std::vector<TObject*> delete_list; 
  };

}


#endif
