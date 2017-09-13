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
  
  const double defaultLevelThresholds[] = {0,1e-10,1e-7,1e-5,1e-3,1e-2,0.05}; //these are CDF thresholds 
  static StereographicGrid defaultSegmentationScheme; 
  static UCorrelator::ConstantPointingResolutionModel defaultPointingResolutionModel; 

  class ProbabilityMap 
  {
    public: 

      /** this encodes configuration for the ProbabilityMap */ 
      struct Params
      {

        Params() 
        {
          seg = &defaultSegmentationScheme; //default
          point = &defaultPointingResolutionModel; 
          n_level_thresholds = sizeof(defaultLevelThresholds) / sizeof(*defaultLevelThresholds); 
          level_cdf_thresholds = defaultLevelThresholds; 
          density_cutoff = 1e-10; 
          projection = BACKWARD; 
          collision_detection = true; 
          dataset = RampdemReader::rampdem; 
        }

        const AntarcticSegmentationScheme * seg; 
        const PointingResolutionModel * point;  
        int n_level_thresholds;
        const double * level_cdf_thresholds; 
        double density_cutoff; 
        RampdemReader::dataSet dataset; 

        enum ProjectionMode
        {
          MC, 
          BACKWARD
        } projection; 

        bool collision_detection; 

        struct CollisionDetectionParams
        {
          CollisionDetectionParams() 
          {
            dx = 1000; 
            grace = 0;
          } 
          double dx; 
          double grace; 
        } collision_params;
       
        struct MCParams
        {
          MCParams() 
          {
            n = 100e6; 
          }
         long long n;

        } mc_params; 

        struct BackwardParams
        {
          BackwardParams() 
          {
            num_samples_per_bin = 64; 
            el_cutoff = 0; 
            random_samples = false; 
          }

          int num_samples_per_bin; 
          double el_cutoff; 
          bool random_samples;
        } backwards_params; 
      }; 

      /** Initialize a probability map. If you pass 0 for params, the defaults are used; 
       */
      ProbabilityMap( const Params  * p  = 0 );


      /** Convert between pdf p-value and cdf p value */ 
      static double cdf2density(double p); 
      static double density2cdf(double density); 



      virtual ~ProbabilityMap() { ; } 

      /** Used to combine information from various maps */ 
      int combineWith(const ProbabilityMap & other); 

      /** Add a point to the clustering. Returns the number of segments that had a non-zero contributions */ 
      int add(const AnitaEventSummary * sum , const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak = 0, double weight = 1, TFile * debugfile = 0); 


      /** This method actually does most of the hard work. 
       *  
       *
       * Some day I'll document it. 
       *
       *  
       *
       */
      void computeContributions(const AnitaEventSummary * sum, const Adu5Pat * pat, 
                                AnitaPol::AnitaPol_t pol, int peak, 
                                std::vector<std::pair<int,double> > & contribution, 
                                std::vector<std::pair<int,double> > * base_contributions = 0,
                                std::vector<std::pair<int,double> > * occlusion = 0, 
                                TFile * debugfile = 0) const; 

      /** Check the overlap of a point with the probability map.
       *  If the point is already in the probability map, remove_self_contribution should be true so that it won't count against itself. If you are checking a point not in the map already with the map, you should set it to false. 
       */ 
      double overlap(const AnitaEventSummary * sum , const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak = 0, 
                     std::vector<std::pair<int,double> >  * overlapped_bases = 0, bool remove_self_contribution =true) const ; 

      /* These are probability densities */ 
      const double* getDensitySums() const { return &ps[0]; } 
      const double* getOccludedFractionSum() const { return &fraction_occluded[0]; } 

      const int* getNAboveLevel(int level) const { return & NAboveLevel[level][0]; } 

      size_t NLevels() const { return p.n_level_thresholds; } 
      double getLevel(int level) const { return p.level_cdf_thresholds[level]; } 

      const AntarcticSegmentationScheme * segmentationScheme() const { return p.seg; } 
      double minDensity() const { return p.density_cutoff; } 

      /* Return the number of bases considered ...  this should match BaseList::getBases() + BaseList::getPaths(
       *  for the right ANITA version
       *
       * The bases are indexed with stationary bases first followed by paths. 
       **/ 
      size_t getNBases() const { return base_ps.size(); } 

      const double  * getBaseDensitySums()  const { return &base_ps[0]; } 
      const int* getBaseNAboveLevel(int level) const { return & baseNAboveLevel[level][0]; } 

      
    private:
      Params p; 

      //indexed by segment
      std::vector<double> ps; 
      std::vector<double> fraction_occluded; 
      //indexed by level then segment 
      std::vector< std::vector<int> > NAboveLevel; 
      std::vector<double> levels_p; 

      //mapping of segment to base in segment  (stationary bases only) 
      std::vector<std::vector<int> > basesInSegment; 

      //indexed by base
      std::vector< std::vector<int> > baseNAboveLevel; 
      std::vector<double> base_ps; 

      ClassDefNV(ProbabilityMap, 3); 
  }; 
}

#endif
