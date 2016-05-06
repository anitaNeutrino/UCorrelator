#ifndef UCORRELATOR_ANALYZER_H
#define UCORRELATOR_ANALYZER_H

#include "Correlator.h" 
#include "WaveformCombiner.h"
#include <vector>
#include "TH2D.h"

#include "AnitaEventSummary.h" 
class FilteredAnitaEvent; 
class UsefulAdu5Pat; 
class AnalysisWaveform; 

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
      const Correlator * getCorrelator() const { return &corr; } 

      /** Retrieve the internal waveform combiner. Note that if multiple polarizations are analyzed, the WaveformCombiner's internal state
       *  will be related to the last polarization used. */ 
      const WaveformCombiner * getWaveformCombiner() const { return &wfcomb; } 

      /** Return the correlation map for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const TH2 * getCorrelationMap(AnitaPol::AnitaPol_t pol) const  { return correlation_maps[pol] ; } 

      /** Return the ith zoomed correlation map for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const TH2 * getZoomedCorrelationMap(AnitaPol::AnitaPol_t pol, int i) const { return zoomed_correlation_maps[pol][i]; }

      /** Return the ith coherent waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getCoherent(AnitaPol::AnitaPol_t pol, int i) const { return coherent[pol][i]; } 

       /** Return the ith deconvolved waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getDeconvolved(AnitaPol::AnitaPol_t pol, int i) const { return deconvolved[pol][i]; } 


    private:

      void fillWaveformInfo(const AnalysisWaveform * wf, AnitaEventSummary::WaveformInfo * info); 
      void fillPointingInfo(double rough_phi, double rough_theta, AnitaEventSummary::PointingHypothesis * point, UsefulAdu5Pat * pat); 
      void fillFlags(const FilteredAnitaEvent * fae, AnitaEventSummary::EventFlags * flags, UsefulAdu5Pat * pat); 

      TH2D* correlation_maps[2]; 
      std::vector<TH2D*>  zoomed_correlation_maps[2]; 
      std::vector<AnalysisWaveform *> coherent[2]; 
      std::vector<AnalysisWaveform *> deconvolved[2]; 

      const AnalysisConfig * cfg; 
      Correlator corr; 
      WaveformCombiner wfcomb; 
      TH2D zoomed; 
      double maprms; 
      FFTtools::DigitalFilter *  power_filter; 
      bool interactive; 
  };

}


#endif
