#ifndef UCORRELATOR_PROBABILITY_MAP
#define UCORRELATOR_PROBABILITY_MAP

/** Probability Map Clustering implementation
 * 
 * For each AntarcticSegment, defined by the AntarcticSegmentationScheme,  we store: 
 *
 *   p: the sum of p-values for a given bin 
 *   NAbove: the number of events with a p-value greater than X in a bin (where X is defined by the level thresholds) 
 *
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu>
 *
 */ 

#include "AntarcticaGeometry.h" 
#include "AnitaEventSummary.h" 
#include "PointingResolutionModel.h" 
#include <vector> 

class Adu5Pat; 
class TFile; 

namespace UCorrelator
{
  class PointingResolutionModel; 
  const double defaultLevelThresholds[] = {0.05, 0.01, 0.001, 0.00001, 0 }; 


  class ProbabilityMap 
  {
    public: 
      /** Initialize a probability map
       * @param seg The segmentation scheme or NULL to use the default (stereographic with defaults) 
       * @param p   The pointing resolution model or NULL to use the default (default constant pointing model) 
       *
       */
      ProbabilityMap(const AntarcticSegmentationScheme * seg= NULL, 
                     const PointingResolutionModel * p = NULL, 
                     int NlevelThresholds = sizeof(defaultLevelThresholds)/sizeof(*defaultLevelThresholds) ,
                     const double * level_thresholds = defaultLevelThresholds ,
                     double prob_cutoff = 1e-6 , 
                     int num_samples_per_bin = 64
                     );

      virtual ~ProbabilityMap() { ; } 

      /** Used to combine information from various maps */ 
      int combineWith(const ProbabilityMap & other); 

      /** Add a point to the clustering. Returns the number of segments that had a non-zero contributions */ 
      int add(const AnitaEventSummary * sum , const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak = 0, double weight = 1, TFile * debugfile = 0); 


      /** This method actually does most of the hard work
       */
      void computeContributions(const AnitaEventSummary * sum, const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak,  std::vector<std::pair<int,double> > & contribution, TFile * debugfile = 0) const; 
      /** Check the overlap of a point with the probability map.
       *  If the point is already in the probability map, remove_self_contribution should be true so that it won't count against itself. If you are checking a point not in the map already with the map, you should set it to false. 
       */ 
      double overlap(const AnitaEventSummary * sum , const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak = 0, bool remove_self_contribution =true) const ; 

      const double* getProbabilities() const { return &ps[0]; } 

      const int* getNAboveLevel(int level) const { return & NAboveLevel[level][0]; } 

      size_t NLevels() const { return levels.size(); } 
      double getLevel(int level) const { return levels[level]; } 

      const AntarcticSegmentationScheme * segmentationScheme() const { return g; } 
      double probabilityCutoff() const { return min_p; } 

    private:
      const AntarcticSegmentationScheme * g; 
      const PointingResolutionModel * prm; 
      std::vector<double> ps; 
      std::vector< std::vector<int>> NAboveLevel; 
      std::vector<double> levels; 
      double min_p; 
      int nsamples; 

      ClassDefNV(ProbabilityMap, 1); 
  }; 
}

#endif
