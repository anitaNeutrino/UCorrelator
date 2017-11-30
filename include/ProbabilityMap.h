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
#include "TMutex.h"

class Adu5Pat; 
class TFile; 

namespace UCorrelator
{
  
  const double defaultLevelThresholds[] = {1,2,3,4,5,6,7,8,9,10}; //these are give in in mahalanobis distance 
  const AntarcticSegmentationScheme &  defaultSegmentationScheme(); 
  const UCorrelator::PointingResolutionModel & defaultPointingResolutionModel(); 

  class ProbabilityMap   : public TObject 
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
          : seg(&defaultSegmentationScheme()), 
            point(&defaultPointingResolutionModel()), 
            level_thresholds(defaultLevelThresholds, defaultLevelThresholds + sizeof(defaultLevelThresholds) / sizeof(*defaultLevelThresholds)), 
            dataset(RampdemReader::rampdem), 
            refract(0), 
            maximum_distance(20),
            min_p_on_continent (1e-3) , 
            projection(BACKWARD), 
            collision_detection(true) , 
            max_dphi(5), 
            max_dtheta(5),
            verbosity(0) 
        {

        }

        const AntarcticSegmentationScheme * seg; 
        const PointingResolutionModel * point;  
        std::vector<double> level_thresholds; 
        RampdemReader::dataSet dataset; 
        const Refraction::Model * refract; 
        double maximum_distance; 
        double min_p_on_continent; 

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
            num_samples_per_bin = 9; 
            enhance_threshold = 0.01; //threshold to enhance if any delauanay region has an area greater than this 
            max_enhance = 4; 
            el_cutoff = 0; 
            random_samples = false; 
          }

          int num_samples_per_bin; 
          double el_cutoff; 
          double enhance_threshold;
          int max_enhance; 
          bool random_samples;
        } backwards_params; 

        double max_dphi; 
        double max_dtheta; 
        int verbosity; 

        ClassDef(Params,6); 
        virtual ~Params() { ; }  
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



      enum OverlapMode
      {
        OVERLAP_SUM_SQRTS,  //use sum of sqrts  
        OVERLAP_SQRT_SUMS,  //use sqrt of sums
        OVERLAP_SUMS,       //use sums (doesn't really make sense) 
        OVERLAP_MAXES       //use maxes (emulating the old clustering style) 
      }; 

      /** Check the overlap of a point with the probability map.
       *  If the point is already in the probability map, remove_self_contribution should be true so that it won't count against itself. If you are checking a point not in the map already with the map, you should set it to false. 
       */ 

      double overlap(const AnitaEventSummary * sum , const Adu5Pat * pat, 
                     AnitaPol::AnitaPol_t pol, int peak = 0, bool normalized = false, 
                     double weight = 1, 
                     std::vector<std::pair<int,double> >  * bases = 0,
                     OverlapMode mode = OVERLAP_SUM_SQRTS, 
                     bool remove_self_contribution =true) const ; 

      /* These are probability sums */ 
      const double* getProbSums(bool normalized = false) const { return normalized ? &ps_norm[0] : &ps[0]; } 
      const double* getProbSumsWithoutBases(int base_level, bool normalized = false) const { return normalized ? &ps_norm_without_base[base_level][0] : &ps_without_base[base_level][0]; } 
      const double* getProbSqrtSums(bool normalized = false) const { return normalized ? &sqrt_ps_norm[0] : &sqrt_ps[0]; } 
      const double* getProbSqrtSumsWithoutBases(int base_level, bool normalized = false) const { return normalized ? &sqrt_ps_norm_without_base[base_level][0] : &sqrt_ps_without_base[base_level][0]; } 
      const double* getProbMaxes(bool normalized = false) const { return normalized ? &max1_ps_norm[0] : & max1_ps[0]; } 
      const double* getProbSecondMaxes(bool normalized = false) const { return normalized ? &max2_ps_norm[0] : & max2_ps[0]; } 
      const double* getOccludedFractionSum() const { return &fraction_occluded[0]; } 

      const int* getNAboveLevel(int level, bool normalized = false) const { return normalized ? &n_above_level_norm[level][0] : &n_above_level[level][0]; } 
      const int* getNAboveLevelWithoutBases(int level, bool normalized = false) const { return normalized ? &n_above_level_without_base_norm[level][0] : &n_above_level_without_base[level][0]; } 

      const double* getWgtAboveLevel(int level, bool normalized=false) const
      { return normalized ? &wgt_above_level_norm[level][0] : & wgt_above_level[level][0]; } 
      const double* getWgtAboveLevelWithoutBases(int level, bool normalized=false) const 
      { return normalized ? &wgt_above_level_without_base_norm[level][0] : &wgt_above_level_without_base[level][0]; } 

      size_t NLevels() const { return p.level_thresholds.size(); } 
      double getLevel(int level) const { return p.level_thresholds[level]; } 

      const AntarcticSegmentationScheme * segmentationScheme() const { return p.seg; } 
      double maxDistance() const { return p.maximum_distance; } 

      /* Return the number of bases considered ...  this should match BaseList::getBases() + BaseList::getPaths(
       *  for the right ANITA version
       *
       * The bases are indexed with stationary bases first followed by paths. 
       **/ 
      size_t getNBases() const { return base_sums.size(); } 

      const double  * getBaseSums(bool normalized = false)  const 
      { return normalized ? &base_sums_norm[0] : &base_sums[0]; } 
      const int* getBaseNAboveLevel(int level, bool normalized = false) const
      { return normalized ? &base_n_above_level_norm[level][0] : & base_n_above_level[level][0]; } 


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

      int makeMultiplicityTable(int level, bool blind = true) const; 
      
    private:
      Params p; 

      //indexed by segment
      std::vector<double> ps; 
      std::vector< std::vector<double> >  ps_without_base; //like ps, but require that no base above level is contained
      std::vector<double> ps_norm; //like ps, but normalized so integral is 1 
      std::vector< std::vector<double> >ps_norm_without_base; //like ps_without_base, but normalized so integral is 1 

      std::vector<double> max1_ps;  //maximum value
      std::vector<double> max1_ps_norm; //max1_ps, but normalized so integral is 1 
      std::vector<double> max2_ps;  //second to maximum value 
      std::vector<double> max2_ps_norm; //max2_ps, but normalized so integral is 1 


      // these are all the sums of the square roots instead, needed for computing the overlaps properly
      std::vector<double> sqrt_ps; 
      std::vector< std::vector<double> >  sqrt_ps_without_base; 
      std::vector<double> sqrt_ps_norm; 
      std::vector< std::vector<double> >sqrt_ps_norm_without_base;
 
      std::vector<double> fraction_occluded; 
      //indexed by level then segment 
      std::vector< std::vector<int> > n_above_level; 
      std::vector< std::vector<int> > n_above_level_norm; 
      std::vector< std::vector<double> > wgt_above_level; // like n_above_level, but 1/N_over_level_per_event is put in each bin, so that the number of contributing events can be reliably determined 
      std::vector< std::vector<double> > wgt_above_level_norm; // like n_above_level, but 1/N_over_level_per_event is put in each bin, so that the number of contributing events can be reliably determined 


      std::vector< std::vector<int> > n_above_level_without_base; 
      std::vector< std::vector<int> > n_above_level_without_base_norm; 
      std::vector< std::vector<double> > wgt_above_level_without_base; 
      std::vector< std::vector<double> > wgt_above_level_without_base_norm; 

      //indexed by base
      std::vector< std::vector<int> > base_n_above_level; ; 
      std::vector< std::vector<int> > base_n_above_level_norm; ; 
      std::vector<double> base_sums; 
      std::vector<double> base_sums_norm; 

      //guards the add method (everything else doesn't touch the internals) 
      TMutex m; 


      ClassDefNV(ProbabilityMap, 11); 
  }; 
}


#endif
