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
 #include "TCanvas.h"

class Adu5Pat; 
class TFile; 

namespace UCorrelator
{
  
  const AntarcticSegmentationScheme &  defaultSegmentationScheme(); 
  const UCorrelator::PointingResolutionModel & defaultPointingResolutionModel(); 

  class ProbabilityMap   : public TObject 
  {
    public: 

      /** Conversions between integrated probability and mahalanobis_distance */
      static double cdf2dist(double p) { return sqrt(-2 * log1p(-p)); }
      static double dist2cdf(double d) { return -expm1(-d*d/2.);      }  

      static double get_two_pi_sqrt_det(double sigma1, double sigma2, double corr) { return  2*M_PI * sigma1 * sigma2 * sqrt(1.-corr*corr); }
      static double get_inv_two_pi_sqrt_det(double sigma1, double sigma2, double corr) { return pow ( get_two_pi_sqrt_det(sigma1,sigma2,corr),  double(-1)); }

      /** Conversion between probability density and mahalanobis  distance */ 
      static double dist2dens(double dist, double inv_two_pi_sqrt_det)  { return inv_two_pi_sqrt_det * exp(-dist*dist/2.); }
      static double dens2dist(double dens, double two_pi_sqrt_det)  {  return  sqrt(-2*log (two_pi_sqrt_det * dens )); }  


      /** this encodes configuration for the ProbabilityMap */ 
      struct Params
      {

        Params() 
          : seg(&defaultSegmentationScheme()), 
            point(&defaultPointingResolutionModel()), 
            dataset(RampdemReader::rampdem), 
            refract(0), 
            maximum_distance(3),// how many sigma to plot on prob map
            // min_p_on_continent (1e-3) , //prob sum of an event on ground should larger than this threshold
            projection(BACKWARD), 
            collision_detection(false) , 
            max_dphi(5), 
            max_dtheta(5),
            radius(20000),
            verbosity(0) 
        {

        }

        const AntarcticSegmentationScheme * seg; 
        const PointingResolutionModel * point;  
        RampdemReader::dataSet dataset; 
        const Refraction::Model * refract; 
        double maximum_distance; 
        double min_p_on_continent;
        double radius; 

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
       /** Used to remove information from various maps */ 
      int removeWith(const ProbabilityMap & other); 
      void maskingWithMap(double & nMasked, double & nNotMasked, const ProbabilityMap & other);

      /** Add a point to the clustering. Returns the number of segments that had a non-zero contributions */ 
      int add(int & NOverlapedBases, double & p_ground, const AnitaEventSummary * sum , const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak = 0, double weight = 1, TFile * debugfile = 0); 



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
      }; 

      /** Check the overlap of a point with the probability map.
       *  If the point is already in the probability map, remove_self_contribution should be true so that it won't count against itself. If you are checking a point not in the map already with the map, you should set it to false. 
       */ 

      double overlap(const AnitaEventSummary * sum , const Adu5Pat * pat, 
                     AnitaPol::AnitaPol_t pol, int peak = 0, bool normalized = false, 
                     double weight = 1, 
                     std::vector<std::pair<int,double> >  * bases = 0,
                     OverlapMode mode = OVERLAP_SUM_SQRTS, 
                     bool remove_self_contribution =true,
                     std::vector<std::pair<int,double> > * segments = 0, //way 
                     std::vector<std::pair<int,double> > * max_dens = 0,  //too 
                     double * inv_two_pi_sqrt_det = 0 //much
                     ) const ;   //auxiliary stuff 

      /* These are probability sums */ 
      double getProbSumsIntegral(bool normalizd = false) const; 
      const double* getProbSums(bool normalized = false) const { 
          std::vector <double> temp_map(p.seg->NSegments(),0);
          temp_map = normalized ? ps_norm : ps; 
          if(blind){
            for (int i =0; i< p.seg->NSegments(); i++){
              if(round(mapOfClusterSizes[i]) == 1 and uniform_ps_without_base[i] != 0){
                temp_map[i] = 0;
              }
            }
            return &temp_map[0];
          }
      } 
      const double* getProbSqrtSums(bool normalized = false) const { return normalized ? &sqrt_ps_norm[0] : &sqrt_ps[0]; } 
      const double* getOccludedFractionSum() const { return &fraction_occluded[0]; } 

      const double* getUniformPS() const
      { return  &uniform_ps[0] ; } 
      const double* getBaseWeightedUniformPS() const 
      { return  &uniform_ps_weighted_by_base[0]; }
      const double* getUniformPSwithoutBase() const 
      { return  &uniform_ps_without_base[0]; } 
      const double* getUniformPSwithBase() const 
      { return  &uniform_ps_with_base[0]; } 
      const double* getClusterSizes() const {
          std::vector <double> temp_map(p.seg->NSegments(),0);
          temp_map = mapOfClusterSizes; 
          const int hpolSignalClusterIndex[] = {3, 7,12,13,34,47,48,29,51,52,50,44,20, 9,24,28,27,38,19,54,46,42,39,41};
          if(blind){
            for (int i =0; i< p.seg->NSegments(); i++){
              // temp_map[i] = 4.0;
              if(round(mapOfClusterSizes[i]) == 1 and uniform_ps_without_base[i] != 0){
                // check whether this singlets is exsit in Hpol cluster indexes
                if ( std::find(std::begin(hpolSignalClusterIndex), std::end(hpolSignalClusterIndex), mapOfClusterIndexs[i]) != std::end(hpolSignalClusterIndex)){
                  temp_map[i] = 1;
                }else{
                  temp_map[i] = 0;
                }
              }else if(round(mapOfClusterSizes[i]) == 1 ){
                temp_map[i] = 1;
              }else if(round(mapOfClusterSizes[i]) > 1 and round(mapOfClusterSizes[i]) < 6){
                temp_map[i] = 2;
              }else if(round(mapOfClusterSizes[i]) >= 6 ){
                temp_map[i] = 3;
              }else{
                temp_map[i] = 0;
              }
            }
            return &temp_map[0];
          }
      }


      const AntarcticSegmentationScheme * segmentationScheme() const { return p.seg; } 
      double maxDistance() const { return p.maximum_distance; } 

      /* Return the number of bases considered ...  this should match BaseList::getBases() + BaseList::getPaths(
       *  for the right ANITA version
       *
       * The bases are indexed with stationary bases first followed by paths. 
       **/ 
      size_t getNBases() const { return eventCountPerBase.size(); } 
      const int* getEventCountPerBase() const
      { return &eventCountPerBase[0]; } 


      /** This will take a set of of values (indexed by segment) 
       *  and form groupings where all adjacent non-zero values are grouped together and the value
       *  is at all of those segments is the integral of the group (stored in counts). Optionally,
       *  will also fill a vector with the distribution of integrals. 
       *
       *  Returns the number of groupings. 
       *
       *  */ 
      int countBasesInThisSegment(int seg) const;
      // int groupAdjacent(const double * vals_to_group, std::vector<std::vector<int> > *  groups  = 0, double* counts = 0, std::vector<double>  * distribution = 0, double val_threshold = 0) const; 
      //return the clusters results to clusterSizes and mapOfClusterSizes
      //threshold is the min p on a bin to do clustering. Set to 0 so it is trival.
      int doClustering(double threshold = 0);
      // given event's longitude and latitude, return the cluster id that this event belongs to.
      //must call after the doClustering
      void evaluateEvent(double & indexOfCluster, double & sizeOfCluster, double theta,double phi,const Adu5Pat * gps);

      int dumpNonZeroBases() const; 
      // use the results form doClustering
      std::pair<int, int> showClusters(int draw = 1, bool blind = true, const char * option = "colz") const; 
      
    private:
      Params p; 

      //indexed by segment
      std::vector<double> ps; 
      std::vector<double> ps_norm; //like ps, but normalized so integral is 1 
      // these are all the sums of the square roots instead, needed for computing the overlaps properly
      std::vector<double> sqrt_ps; 
      std::vector<double> sqrt_ps_norm; 
 
      std::vector<double> fraction_occluded; 

      // using a uniform distribution to put into each segment for each events
      std::vector<double> uniform_ps; 
      // if an event does not overlap with anybase, then the uniform distribution sum to 0.000001 instead of 1.
      //Since it is non-zero, all events will be clustered but when we sum the uniform_ps, only the events near base will be count. 
      std::vector<double> uniform_ps_weighted_by_base; 
      std::vector<double> uniform_ps_with_base; 
      std::vector<double> uniform_ps_without_base; 

      //indexed of Base, number of events near this base.
      std::vector<int> eventCountPerBase;

      std::vector<double> mapOfClusterSizes;
      std::vector<double> mapOfClusterIndexs;
      std::vector<double> clusterSizes;
      std::vector<double> clusterSizes_with_base;
      std::vector<double> clusterSizes_without_base;
      bool blind = false;

      //guards the add method (everything else doesn't touch the internals) 
      TMutex m; 


      ClassDefNV(ProbabilityMap, 12); 
  }; 
}


#endif
