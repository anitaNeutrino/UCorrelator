#ifndef UCORRELATOR_GUI
#define UCORRELATOR_GUI

/** Various GUI tools 
 *
 */

#include "TH2.h" 
#include "TMultiGraph.h" 
#include "TPaveText.h" 
#include "TCanvas.h" 
#include "TGaxis.h" 
#include <vector>
#include "PeakFinder.h" 
#include "TMarker.h" 
#include "TEllipse.h" 

class AnitaEventSummary ; 
class FilteredAnitaEvent; 
class AnalysisWaveform; 

namespace UCorrelator
{
  class WaveformCombiner; 
  class Analyzer; 

  namespace gui
  {
    
    class Map : public TH2D 
    {
      public:
        /* Construct the map.
         *
         * The Summary doesn't have to be filled yet except for the pulsers / MC truth */ 

        Map(const TH2D & hist,const FilteredAnitaEvent* f,  WaveformCombiner * comb, WaveformCombiner * comb_filtered, AnitaPol::AnitaPol_t pol, const AnitaEventSummary * sum); 
        Map(const Map & other); 
        virtual ~Map(); 
        virtual void Paint(Option_t * opt = ""); 
        virtual void ExecuteEvent(int event, int x, int y); 
        void SetUseFiltered(); // *MENU*
        void SetUseUnfiltered(); // *MENU*
        void addRough(double x, double y); 
        void addRough(const std::vector<std::pair<double,double> > & rough); 
        void closeCanvas(); 
        void drawWf(double theta,double phi); 
        void addFine(const AnitaEventSummary::PointingHypothesis & p); 
        void clear(); 

      private:
        std::vector<TMarker> rough_m; 
        std::vector<TMarker> fine_m; 
        std::vector<TEllipse> fine_e; 
        std::vector<TMarker> specials; 
        TCanvas * wfpad; 
        const FilteredAnitaEvent *f; 
        WaveformCombiner *c; 
        WaveformCombiner *cf; 
        TMarker * clicked; 
        bool use_filtered; 
        AnitaPol::AnitaPol_t pol;
        AnalysisWaveform * coherent; 
        AnalysisWaveform * deconvolved; 
        double last_theta, last_phi; 
        TGaxis  heading_axis; 

        ClassDef(Map,2); 
    }; 


    /* TODO clickable stuff
    class Waveform: public TMultiGraph  
    {



    }; */


    class SummaryText : public TPaveText
    { 
      public: 
        SummaryText(int i, AnitaPol::AnitaPol_t pol, const Analyzer* a, int use_filtered); 

        ClassDef(SummaryText,1); 
    }; 

  }






}

#endif
