#ifndef UCORRELATOR_PROBABILITY_MAP
#define UCORRELATOR_PROBABILITY_MAP

/** Probability Map Clustering implementation
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu>
 *
 */ 

#include "AntarcticaGeometry.h" 
#include "AnitaEventSummary.h" 
#include "PointingResolutionModel.h" 

namespace UCorrelator
{

  const double defaultLevelThresholds[] = {0.05, 0.01, 0.001, 0.00001, 0 }; 


  class ProbabilityMap 
  {
    public: 
      ProbabilityMap(const AntarcticSegmentationScheme * seg, 
                     int NlevelThresholds = sizeof(defaultLevelThresholds)/sizeof(*defaultLevelThresholds) ,
                     const double * level_thresholds = defaultLevelThresholds ); 

      void combineWith(const ProbabilityMap & other); 
      void add(const AnitaEventSummary * sum , AnitaPol::AnitaPol_t pol, int peak = 0); 
      double  overlap(const AnitaEventSummary * sum , AnitaPol::AnitaPol_t pol, int peak = 0); 


    private:
      const AntarcticSegmentationScheme * g; 
      std::vector<double> ps; 
      std::vector< std::vector<double> > NAboveLevel; 
      std::vector<double> levels; 
  }; 
}

#endif
