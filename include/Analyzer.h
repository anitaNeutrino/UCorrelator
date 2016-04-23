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

  class Analyzer
  {

    public:
      Analyzer(const AnalysisConfig * cfg = 0, bool interactive = false); 
      virtual ~Analyzer();



      void analyze(const FilteredAnitaEvent * fae, AnitaEventSummary *summary); 

//      AnitaEventSummary * analyze(AnitaDataset * dataset, FilterStrategy * strategy, AnitaEventSummary * summary = 0, int eventNumber = -1, int entryNumber=-1); 

      const Correlator * getCorrelator() const { return &corr; } 
      const WaveformCombiner * getWaveformCombiner() const { return &wfcomb; } 

      const TH2 * getCorrelationMap(AnitaPol::AnitaPol_t pol) const  { return correlation_maps[pol] ; } 
      const TH2 * getZoomedCorrelationMap(AnitaPol::AnitaPol_t pol, int i) const { return zoomed_correlation_maps[pol][i]; }
      const AnalysisWaveform * getCoherent(AnitaPol::AnitaPol_t pol, int i) const { return coherent[pol][i]; } 
      const AnalysisWaveform * getDeconvolved(AnitaPol::AnitaPol_t pol, int i) const { return deconvolved[pol][i]; } 



    private:

      void fillWaveformInfo(const AnalysisWaveform * wf, AnitaEventSummary::WaveformInfo * info); 
      void fillPointingInfo(double rough_phi, double rough_theta, AnitaEventSummary::PointingHypothesis * point, UsefulAdu5Pat * pat); 
      void fillFlags(const FilteredAnitaEvent * fae, AnitaEventSummary::EventFlags * flags); 

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
