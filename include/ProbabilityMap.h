#ifndef UCORRELATOR_PROBABILITY_MAP
#define UCORRELATOR_PROBABILITY_MAP

/** Probability Map Clustering implementation
 * 
 * For each AntarcticSegment, defined by the AntarcticSegmentationScheme,  we store: 
 *
 *
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu>
 *
 */ 

#include <math.h> 
#include "AntarcticaGeometry.h" 
#include "AnitaEventSummary.h" 
#include "RefractionModel.h" 
#include "PointingResolutionModel.h" 
#include <vector> 
#include "TBits.h" 

class Adu5Pat; 
class TFile; 

namespace UCorrelator
{
  
  const double defaultLevelThresholds[] = {1,2,3,4,5,6,7,8,9,10}; //these are give in in mahalanobis distance 
  const AntarcticSegmentationScheme &  defaultSegmentationScheme(); 
  const UCorrelator::PointingResolutionModel & defaultPointingResolutionModel(); 

  class ProbabilityMap  
  {
    public: 

      /** Conversions between integrated probability and mahalanobis_distance */
      static double cdf2dist(double p) { return sqrt(-2 * log1p(-p)); }
      static double dist2cdf(double d) { return -expm1(-d*d/2.);      }  

      static double get_two_pi_sqrt_det(double sigma1, double sigma2, double corr) { return  2*M_PI * sigma1 * sigma2 * sqrt(1.-corr*corr); }
      static double get_inv_two_pi_sqrt_det(double sigma1, double sigma2, double corr) { return pow ( get_two_pi_sqrt_det(sigma1,sigma2,corr),   -1); }

      /** Conversion between probability density and mahalanobis  distance */ 
      static double dist2dens(double dist, double inv_two_pi_sqrt_det)  { return inv_two_pi_sqrt_det * exp(-dist*dist/2.); }
      static double dens2dist(double dens, double two_pi_sqrt_det)  {  return  sqrt(-2*log (two_pi_sqrt_det * dens )); }  


      /** this encodes configuration for the ProbabilityMap */ 
      struct Params
      {

        Params() 
        {
          seg = &defaultSegmentationScheme(); //default
          point = &defaultPointingResolutionModel(); 
          n_level_thresholds = sizeof(defaultLevelThresholds) / sizeof(*defaultLevelThresholds); 
          level_thresholds = defaultLevelThresholds; 
          maximum_distance = 10; 
          projection = BACKWARD; 
          collision_detection = true; 
          refract = 0; 
          dataset = RampdemReader::rampdem; 
        }

        const AntarcticSegmentationScheme * seg; 
        const PointingResolutionModel * point;  
        int n_level_thresholds;
        const double * level_thresholds; //[n_level_thresholds] 
        RampdemReader::dataSet dataset; 
        const Refraction::Model * refract; 
        double maximum_distance; 

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
      double computeContributions(const AnitaEventSummary * sum, const Adu5Pat * pat, 
                                AnitaPol::AnitaPol_t pol, int peak, 
                                std::vector<std::pair<int,double> > & contribution, 
                                std::vector<std::pair<int,double> > * base_contributions = 0,
                                std::vector<std::pair<int,double> > * occlusion = 0, 
                                std::vector<std::pair<int,double> > * maximum_density = 0, 
                                TFile * debugfile = 0) const; 

      /** Check the overlap of a point with the probability map.
       *  If the point is already in the probability map, remove_self_contribution should be true so that it won't count against itself. If you are checking a point not in the map already with the map, you should set it to false. 
       */ 
      double overlap(const AnitaEventSummary * sum , const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak = 0, 
                     std::vector<std::pair<int,double> >  * overlapped_bases = 0, bool remove_self_contribution =true) const ; 

      /* These are probability densities */ 
      const double* getDensitySums() const { return &ps[0]; } 
      const double* getDensitySumsWithBases() const { return &ps_with_base[0]; } 
      const double* getDensitySumsWithoutBases() const { return &ps_without_base[0]; } 

      const double* getDensitySumsNormalized() const { return &ps_norm[0]; } 
      const double* getDensitySumsNormalizedWithBases() const { return &ps_norm_with_base[0]; } 
      const double* getDensitySumsNormalizedWithoutBases() const { return &ps_norm_without_base[0]; } 

      const double* getOccludedFractionSum() const { return &fraction_occluded[0]; } 

      const int* getNAboveLevel(int level) const { return & n_above_level[level][0]; } 
      const int* getNAboveLevelWithBases(int level) const { return & n_above_level_with_base[level][0]; } 
      const int* getNAboveLevelWithoutBases(int level) const { return & n_above_level_without_base[level][0]; } 

      const double* getWgtAboveLevel(int level) const { return & wgt_above_level[level][0]; } 
      const double* getWgtAboveLevelWithBases(int level) const { return &wgt_above_level_with_base[level][0]; } 
      const double* getWgtAboveLevelWithoutBases(int level) const { return &wgt_above_level_without_base[level][0]; } 

      size_t NLevels() const { return p.n_level_thresholds; } 
      double getLevel(int level) const { return p.level_thresholds[level]; } 

      const AntarcticSegmentationScheme * segmentationScheme() const { return p.seg; } 
      double maxDistance() const { return p.maximum_distance; } 

      /* Return the number of bases considered ...  this should match BaseList::getBases() + BaseList::getPaths(
       *  for the right ANITA version
       *
       * The bases are indexed with stationary bases first followed by paths. 
       **/ 
      size_t getNBases() const { return base_ps.size(); } 

      const double  * getBaseDensitySums()  const { return & base_ps[0]; } 
      const int* getBaseNAboveLevel(int level) const { return & base_n_above_level[level][0]; } 


      /** This will take a set of of values (indexed by segment) 
       *  and form groupings where all adjacent non-zero values are grouped together and the value
       *  is at all of those segments is the integral of the group (stored in counts). Optionally,
       *  will also fill a vector with the distribution of integrals. 
       *
       *  Returns the number of groupings. 
       *
       *  */ 
      int groupAdjacent(const double * vals_to_group, double* counts, std::vector<double>  * distribution = 0) const; 

      int dumpNonZeroBases() const; 

      
    private:
      Params p; 

      //indexed by segment
      std::vector<double> ps; 
      std::vector<double> ps_with_base; //like ps, but require that a base is contained
      std::vector<double> ps_without_base; //like ps, but require that no base is contained
      std::vector<double> ps_norm; //like ps, but normalized so integral is 1 
      std::vector<double> ps_norm_with_base; //like ps_with_base, but normalized so integral is 1 
      std::vector<double> ps_norm_without_base; //like ps_without_base, but normalized so integral is 1 
 
      std::vector<double> fraction_occluded; 
      //indexed by level then segment 
      std::vector< std::vector<int> > n_above_level; 
      std::vector< std::vector<double> > wgt_above_level; // like n_above_level, but 1/N_over_level_per_event is put in each bin, so that the number of contributing events can be reliably determined 



      //mapping of segment to base in segment  (for stationary bases) 
      std::vector<std::vector<int> > bases_in_segment; 

      //This stores the number of events which both have both at least this level with the segment AND with a base, so can be used as a proxy for if a base is present or not
      std::vector< std::vector<int> > n_above_level_with_base; 
      std::vector< std::vector<double> > wgt_above_level_with_base; 

      std::vector< std::vector<int> > n_above_level_without_base; 
      std::vector< std::vector<double> > wgt_above_level_without_base; 


      //indexed by base
      std::vector< std::vector<int> > base_n_above_level; ; 
      std::vector<double> base_ps; 


      ClassDefNV(ProbabilityMap, 7); 
  }; 
}


#endif
