#ifndef UCORRELATOR_ANALYZER_H
#define UCORRELATOR_ANALYZER_H

#include "Correlator.h" 
#include "WaveformCombiner.h"
#include "UCorrelatorGUI.h" 
#include <vector>
#include "TH2D.h"
#include "simpleStructs.h"
#include "AnitaEventSummary.h" 
#include "ResponseManager.h" 
#include "UsefulAnitaEvent.h" 

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



      /** Analyze the event and put the results into the summary. The summary is reinitialized, clobbering everything.  MCTruth may be passed to initalize truth part of summary. */
      void analyze(const FilteredAnitaEvent * event, AnitaEventSummary *summary, const TruthAnitaEvent * truth = 0); 

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
      const gui::Map * getCorrelationMap(AnitaPol::AnitaPol_t pol) const  { return correlation_maps[pol] ; } 

      /** Return the ith zoomed correlation map for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const gui::Map * getZoomedCorrelationMap(AnitaPol::AnitaPol_t pol, int i) const { return zoomed_correlation_maps[pol][i]; }

      /** Return the ith coherent waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getCoherent(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return coherent[pol][filtered ? 1 : 0][i]; } 

       /** Return the ith deconvolved waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getDeconvolved(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return deconvolved[pol][filtered ? 1: 0][i]; } 

      /** Return the ith coherent xpol waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getCoherentXpol(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return coherent_xpol[pol][filtered ? 1: 0][i]; } 

       /** Return the ith deconvolved xpol waveform for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const AnalysisWaveform * getDeconvolvedXpol(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return deconvolved_xpol[pol][filtered ? 1: 0][i]; } 

      /** Return the ith coherent averaged power for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const TGraphAligned * getCoherentPower(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return coherent_power[pol][filtered ? 1: 0][i]; } 

       /** Return the ith deconvolved averaged power for the polarization. For this to work, Analyzer must have been constructed with interactive= true and the polarization asked for must have been enabled in the config. */ 
      const TGraphAligned * getDeconvolvedPower(AnitaPol::AnitaPol_t pol, int i, bool filtered = false) const { return deconvolved_power[pol][filtered ? 1: 0][i]; } 

      AnitaResponse::ResponseManager * getResponseManager() { return &responses; } 

      /* Return the summary. This only works in interactive mode. */ 
      const AnitaEventSummary* getSummary() const { return &last; }


      /** Populate the pads with a summary of the pointing. Only makes sense if interactive mode is on. If the analyzer wasn't instructed to do 
       * both polarities, then it won't populate any it wasn't instructed to do. If 0 or NULL is passed, a new canvas is made. Added a flag to allow you to draw and populate the info from the filtered waveforms (off by default) */ 
      void drawSummary(TPad *chpol = 0, TPad * cvpol = 0, int draw_filtered = 0) const; 

      double getRoughPhi(AnitaPol::AnitaPol_t pol, int i) const { return rough_peaks[pol][i].first; }
      double getRoughTheta(AnitaPol::AnitaPol_t pol, int i) const { return -rough_peaks[pol][i].second; }

      void clearInteractiveMemory(double frac = 0.5) const; 
      
      /** Sets disallowed antennas with a bitmask */
      void setDisallowedAntennas(ULong64_t hpol=0, ULong64_t vpol=0) {disallowedAnts[0] = hpol; disallowedAnts[1] = vpol; } 
      
      /** Excludes the range for finding peaks in phi (full range is 0 -> 360)*/
      void setExcludePhiRange(double phiMin, double phiMax) {phiRange[0] = phiMin; phiRange[1] = phiMax; exclude = true; }
      /** Excludes the range for finding peaks in theta (full range is -60 -> 40)*/
      void setExcludeThetaRange(double thetaMin, double thetaMax) {thetaRange[0] = thetaMin; thetaRange[1] = thetaMax; exclude = true; }
      /** Excludes the phi and theta range for finding peaks*/
      void setExcludeThetaPhiRange(double thetaMin, double thetaMax, double phiMin, double phiMax) {thetaRange[0] = thetaMin; thetaRange[1] = thetaMax; phiRange[0] = phiMin; phiRange[1] = phiMax; exclude = true; }
      
      /** Looks only at the specified phi range for finding peaks (full phi range is 0 -> 360)*/
      void setPhiRange(double phiMin, double phiMax) {phiRange[0] = phiMin; phiRange[1] = phiMax; exclude = false; }
      /** Looks only at the specified theta range for finding peaks (full theta range is -60 -> 40)*/
      void setThetaRange(double thetaMin, double thetaMax) {thetaRange[0] = thetaMin; thetaRange[1] = thetaMax; exclude = false; }
      /** Looks only at the specified range in phi and theta for finding peaks*/
      void setThetaPhiRange(double thetaMin, double thetaMax, double phiMin, double phiMax) {thetaRange[0] = thetaMin; thetaRange[1] = thetaMax; phiRange[0] = phiMin; phiRange[1] = phiMax; exclude = false; }
      /** Set to look exclusively around a source for peaks or to exclude the angular area around a source from peak finding */
      void setTrackSource(double setLon, double setLat, double setAlt, double setdTheta = 2.5, double setdPhi = 5., bool blockOut = false) {sourceLon = setLon; sourceLat = setLat; sourceAlt = setAlt; dTheta = setdTheta; dPhi = setdPhi; exclude = blockOut; trackSource = true; }
      /** Tracks the sun for peak finding, can be set to exclude or include the sun */
      void setTrackSun(double setdTheta = 2.5, double setdPhi = 5., bool blockOut = false) {dTheta = setdTheta; dPhi = setdPhi; exclude = blockOut; trackSun = true; }
      /** Tracks WAIS for peak finding, include or exclude  */
      void setTrackWAIS(double setdTheta = 2.5, double setdPhi = 5., bool blockOut = false) {sourceLon = AnitaLocations::getWaisLongitude(); sourceLat = AnitaLocations::getWaisLatitude(); sourceAlt = AnitaLocations::getWaisAltitude(); dTheta = setdTheta; dPhi = setdPhi; exclude = blockOut; trackSource = true; }
      /** Tracks LDB for peak finding  */
      void setTrackLDB(double setdTheta = 2.5, double setdPhi = 5., bool blockOut = false) {sourceLon = AnitaLocations::LONGITUDE_LDB; sourceLat = AnitaLocations::LATITUDE_LDB; sourceAlt = AnitaLocations::ALTITUDE_LDB; dTheta = setdTheta; dPhi = setdPhi; exclude = blockOut; trackSource = true; }
			/** Allows you to set extra filters used only for combining waveforms */
			void setExtraFilters(FilterStrategy* extra);
			/** Allows you to set extra filters used only for combining deconvolved waveforms */
			void setExtraFiltersDeconvolved(FilterStrategy* extra);
      
  private:

      void fillWaveformInfo(const AnalysisWaveform * wf, const AnalysisWaveform * xpol_wf, const TGraph * power, AnitaEventSummary::WaveformInfo * info, AnitaPol::AnitaPol_t pol, double rms); 
      void fillPointingInfo(double rough_phi, double rough_theta, AnitaEventSummary::PointingHypothesis * point,
                            UsefulAdu5Pat * pat, double hwAngle, UShort_t triggered_sectors, UShort_t masked_sectors, UShort_t triggered_sectors_xpol, UShort_t masked_sectors_xpol); 
      void fillFlags(const FilteredAnitaEvent * fae, AnitaEventSummary* summary, UsefulAdu5Pat * pat); 
      void fillChannelInfo(const FilteredAnitaEvent* event, AnitaEventSummary* summary);

      gui::Map* correlation_maps[2]; 
      std::vector<gui::Map*>  zoomed_correlation_maps[2]; 
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
      AnitaResponse::ResponseManager responses; 
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
      ULong64_t disallowedAnts[2];
      double phiRange[2];
      double thetaRange[2];
      bool exclude;
      bool trackSource;
      double sourceLon;
      double sourceLat;
      double sourceAlt;
      bool trackSun;
      double dPhi;
      double dTheta;

      mutable std::vector<TObject*> delete_list; 
  };

}


#endif
