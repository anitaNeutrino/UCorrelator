#include "ProbabilityMap.h" 
#include "PointingResolutionModel.h"
#include "Adu5Pat.h" 
#include "assert.h" 
#include "BaseList.h" 
#include "Math/ProbFunc.h"
#include "UsefulAdu5Pat.h" 
 

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
 #define HAVE_DELAUNAY 1
#endif


#ifdef HAVE_DELAUNAY
#include "Math/Delaunay2D.h"
#endif

#include "TGraph2D.h" 
#include "TFile.h" 


ClassImp(UCorrelator::ProbabilityMap); 

const AntarcticSegmentationScheme &  UCorrelator::defaultSegmentationScheme() { static StereographicGrid g; return g; } 
const UCorrelator::PointingResolutionModel & UCorrelator::defaultPointingResolutionModel()  { static UCorrelator::ConstantPointingResolutionModel m; return m; } 

static UCorrelator::ProbabilityMap::Params default_params; 

UCorrelator::ProbabilityMap::ProbabilityMap(const Params * par) 
  :  
    p( par ? *par : default_params), 
    ps(p.seg->NSegments(),0), 
    ps_norm(p.seg->NSegments(),0), 
    sqrt_ps(p.seg->NSegments(),0), 
    sqrt_ps_norm(p.seg->NSegments(),0), 
    fraction_occluded(p.seg->NSegments(), 0), 
    uniform_ps(p.seg->NSegments(),0),
    uniform_ps_weighted_by_base(p.seg->NSegments(),0),
    uniform_ps_with_base(p.seg->NSegments(),0),
    uniform_ps_without_base(p.seg->NSegments(),0),
    mapOfClusterSizes(p.seg->NSegments(),0),
    mapOfClusterIndexs(p.seg->NSegments(),0),
    eventCountPerBase(BaseList::getNumBases() + BaseList::getNumPaths(), 0)
{
  //blank, initialization is above
}


int UCorrelator::ProbabilityMap::add(int& NOverlapedBases, double & p_ground, const AnitaEventSummary * sum, const Adu5Pat * pat, AnitaPol::AnitaPol_t pol,
                                     int peak, double weight, TFile * debugfile) 
{

  std::vector<std::pair<int,double> > segments_to_fill; 
  std::vector<std::pair<int,double> > base_ps_to_fill; 
  std::vector<std::pair<int,double> > occluded_to_fill; 
  std::vector<std::pair<int,double> > max_densities; ; 
  //returned all segments to fill depends on the direction and the max distance.
  //seg_ps_to_fill: return the probability sums for the segments.
  //base_ps_to_fill: return the prob density(not ps) for base.
  //occlude to fill: return the fraction of occluded sample for the segment.
  //max_density_to_fill: return the max prob density for the segment.
  double inv_two_pi_sqrt_det = computeContributions(sum,pat,pol,peak,segments_to_fill, &base_ps_to_fill, &occluded_to_fill, &max_densities, debugfile); 

  if (!inv_two_pi_sqrt_det) return 0; 
  
  TLockGuard lock(&m); 

  NOverlapedBases = base_ps_to_fill.size(); 

  int NOfSegments = segments_to_fill.size(); 
  double norm = 0; 
  for (int i = 0; i < NOfSegments; i++)
  {
    //add all segments' prob sum together, which is the norm factor.
    norm += segments_to_fill[i].second; 
  }
  p_ground = norm;
  // inverse of norm. if normal is too small, let invnorm be 0.
  double invnorm = norm == 0 ? 0 : 1./norm;
  if (p.verbosity > 2) printf("invnorm: %g\n", invnorm); 

  for (int i = 0; i < NOverlapedBases; i++)
  {
    int ibase = base_ps_to_fill[i].first; 
    double dens_base = base_ps_to_fill[i].second * weight; 
    eventCountPerBase[ibase]++; 
  }

  for (int i = 0; i < NOfSegments; i++)
  {
    //loop through each segments that have projection.
    int iseg = segments_to_fill[i].first; 
    double pseg = segments_to_fill[i].second * weight; 
    double pseg_norm = pseg * invnorm * weight; 

    //get the largest and second largest prob density sum for each segment.
    //do it for normalized.

    //prob density sum for this segment. normlized. sqrt. sqrt_normlized
    ps[iseg] += pseg; 
    ps_norm[iseg] += (pseg_norm); 
    sqrt_ps[iseg]  += weight < 0 ? -sqrt(-pseg) : sqrt(pseg); 
    sqrt_ps_norm[iseg] +=  weight < 0 ? -sqrt(-pseg_norm) : sqrt(pseg_norm); 
    //fraction excluded sample. writen to file.
    fraction_occluded[iseg] += occluded_to_fill[i].second; 
    //1/number of segment. Like use uniforms distribution to replace the gaussian. The cut of angular is still using gaussian.
    //no matter a event is near a base or not, it will be added to this map.
    uniform_ps[iseg]+= 1/ double(NOfSegments)* weight;
    if(NOverlapedBases!=0){ 
      uniform_ps_weighted_by_base[iseg]+= 1/ double(NOfSegments)* weight;
      uniform_ps_with_base[iseg]+= 1/ double(NOfSegments)* weight;
    }else{
      uniform_ps_weighted_by_base[iseg]+= 0.000001/ double(NOfSegments)* weight;
      uniform_ps_without_base[iseg]+= 1/ double(NOfSegments)* weight;
    } 
  }
  return NOfSegments; 
}

int UCorrelator::ProbabilityMap::dumpNonZeroBases()  const 
{
  int n = 0;
  for (int i = 0; i < (int) getNBases(); i++)
  {
    if(eventCountPerBase[i] > 0) 
    {
      std::cout <<"BaseIndex = "<< i << "\t BaseName = "<< BaseList::getAbstractBase(i).getName() << "\t eventCount"<<eventCountPerBase[i] << std::endl;; 
      n++; 
    }
  }
 
  return n; 
}

double UCorrelator::ProbabilityMap::overlap(const AnitaEventSummary * sum , const Adu5Pat * pat,  AnitaPol::AnitaPol_t pol, int peak, 
                                            bool normalized, double weight,  std::vector<std::pair<int,double> > * bases,
                                            OverlapMode mode,  bool remove_self, std::vector<std::pair<int, double > > * segments, 
                                            std::vector<std::pair<int,double> > * max_dens, double * inv_two_pi_sqrt_det ) const
{
  weight = 1;
  //so no weight
  std::vector<std::pair<int,double> > segs; 

  double inv = computeContributions(sum,pat,pol,peak,segs,bases,0,max_dens); 
  
  int N = segs.size(); 
  double answer = 0;
  double invnorm = 1;
  if (normalized)
  {
    double norm = 0;
    for (int i =0; i< N; i++)
    {
      norm += segs[i].second; 
    }
    // std::cout<<" norm = " <<norm<<std::endl;
    if (!inv || norm == 0)
    {
      printf("the norm on ground is: %g!\n", norm); 
      return -1; // there can be no overlap 
    }
    invnorm = 1./norm; 
  }

  if (inv_two_pi_sqrt_det) *inv_two_pi_sqrt_det = inv; 


  //pick the right thing to overlap with 
  const double * the_rest = mode == OVERLAP_SUM_SQRTS ? getProbSqrtSums(normalized) : getProbSums(normalized);
 for (int i =0; i< N; i++)
  {
    int seg = segs[i].first; 

    double p_this  =  weight*segs[i].second * invnorm; 

    double p_other = the_rest[seg]; 
    // blild to unknown base singlets to calculate mc overlaping
    if(blind and round(mapOfClusterSizes[i]) == 1 and uniform_ps_without_base[i] != 0){
      p_other = 0; 
    }

    //remove the current event's ps from all segments.
    if (remove_self && mode == OVERLAP_SQRT_SUMS)
    {
      p_other -= p_this; 
    }

    if (remove_self && mode == OVERLAP_SUM_SQRTS) 
    {
      p_other -= sqrt(p_this);  //this might help with floating point precision to do first... maybe. 
    }
    double danswer = (mode == OVERLAP_SUM_SQRTS ? sqrt(p_this) : p_this) * p_other; 
    answer += danswer; 

  }

  if (segments) *segments = segs; 

  // avoid negative roundoff 
  if (answer < 0) answer = 0; 
  return answer; 
}

int UCorrelator::ProbabilityMap::combineWith(const ProbabilityMap & other) 
{

  TString our_scheme; 
  TString their_scheme; 

  segmentationScheme()->asString(&our_scheme); 
  other.segmentationScheme()->asString(&their_scheme); 

  if (our_scheme != their_scheme)
  {
    fprintf(stderr,"Cannot combine ProbabilityMap's with different segmentation schemes (ours: %s, theirs: %s)\n", our_scheme.Data(), their_scheme.Data()); 
    return 1; 
  }
  bool combine_bases = true; 
  if (getNBases() != other.getNBases())
  {
    fprintf(stderr,"Cannot combine bases from ProbabilityMap's since there are different numbers of bases (ours: %d, theirs: %d)\n", (int) getNBases(),(int)  other.getNBases()); 

    combine_bases = false; 
  }


  if (maxDistance() != other.maxDistance())
  {
    fprintf(stderr,"Combining ProbabilityMap's with different max mahalanobis distance (ours=%g,theirs=%g). Continuing, but may not be what you want!\n", maxDistance(), other.maxDistance());  
  }
  //if we made it this far, all is good (except for what hasn't been implemented yet) !
  for (int i = 0; i < segmentationScheme()->NSegments(); i++) 
  {
    ps[i] += other.getProbSums()[i]; 
    ps_norm[i] += other.getProbSums(true)[i]; 
    sqrt_ps[i] += other.getProbSqrtSums()[i]; 
    sqrt_ps_norm[i] += other.getProbSqrtSums(true)[i];
    uniform_ps[i] += other.getUniformPS()[i]; 
    uniform_ps_weighted_by_base[i] += other.getBaseWeightedUniformPS()[i];       
    uniform_ps_with_base[i] += other.getUniformPSwithBase()[i];       
    uniform_ps_without_base[i] += other.getUniformPSwithoutBase()[i];       
  }
  if (combine_bases) 
  {
    for (int i = 0; i < getNBases(); i++)
    {
      eventCountPerBase[i] += other.getEventCountPerBase()[i]; 
    }
  }

  return 0; 
}

int UCorrelator::ProbabilityMap::removeWith(const ProbabilityMap & other) 
{

  TString our_scheme; 
  TString their_scheme; 

  segmentationScheme()->asString(&our_scheme); 
  other.segmentationScheme()->asString(&their_scheme); 

  if (our_scheme != their_scheme)
  {
    fprintf(stderr,"Cannot combine ProbabilityMap's with different segmentation schemes (ours: %s, theirs: %s)\n", our_scheme.Data(), their_scheme.Data()); 
    return 1; 
  }

  bool combine_bases = true; 
  if (getNBases() != other.getNBases())
  {
    fprintf(stderr,"Cannot combine bases from ProbabilityMap's since there are different numbers of bases (ours: %d, theirs: %d)\n", (int) getNBases(),(int)  other.getNBases()); 

    combine_bases = false; 
  }


  if (maxDistance() != other.maxDistance())
  {
    fprintf(stderr,"Combining ProbabilityMap's with different max mahalanobis distance (ours=%g,theirs=%g). Continuing, but may not be what you want!\n", maxDistance(), other.maxDistance());  
  }
  //if we made it this far, all is good (except for what hasn't been implemented yet) !
  for (int i = 0; i < segmentationScheme()->NSegments(); i++) 
  {
    ps[i] -= other.getProbSums()[i]; 
    ps_norm[i] -= other.getProbSums(true)[i]; 
    sqrt_ps[i] -= other.getProbSqrtSums()[i]; 
    sqrt_ps_norm[i] -= other.getProbSqrtSums(true)[i]; 
    uniform_ps[i] -= other.getUniformPS()[i]; 
    uniform_ps_weighted_by_base[i] -= other.getBaseWeightedUniformPS()[i]; 
    uniform_ps_with_base[i] -= other.getUniformPSwithBase()[i];       
    uniform_ps_without_base[i] -= other.getUniformPSwithoutBase()[i];       
  }
  if (combine_bases) 
  {
    for (size_t i = 0; i < getNBases(); i++) 
    {
      eventCountPerBase[i] -= other.getEventCountPerBase()[i]; 
    }
  }
  return 0; 
}
void UCorrelator::ProbabilityMap::maskingWithMap(double & nMasked, double & nNotMasked, const ProbabilityMap & other){
  //if we made it this far, all is good (except for what hasn't been implemented yet) !
  for (int i = 0; i < segmentationScheme()->NSegments(); i++) 
  {
    int hit = 0;
    if(other.getUniformPS()[i]!=0){
      hit = 1;
    }else{
      std::vector<int> new_neighbors;
      segmentationScheme()->getNeighbors(i, &new_neighbors); // find the 8 neighbor seg around this seg. The return vector is in new_neighbors.
      for (size_t j = 0; j < new_neighbors.size(); j++){
        if(other.getUniformPS()[new_neighbors[j]]!=0){ 
          hit = 1;
          break;
        }
      }
    }
    if(hit == 1){
      nMasked += uniform_ps[i];
    }else{
      nNotMasked += uniform_ps[i];
    }
  }

  return ;
}


double  UCorrelator::ProbabilityMap::computeContributions(const AnitaEventSummary * sum, const Adu5Pat * gps, 
                                                       AnitaPol::AnitaPol_t pol, int peak,
                                                      std::vector<std::pair<int,double> > & contribution, 
                                                      std::vector<std::pair<int,double> > * base_contribution, 
                                                      std::vector<std::pair<int,double> > * occlusion, 
                                                      std::vector<std::pair<int,double> > * max_densities, 
                                                      TFile * debugfile) const 
{

  /* compute pointing resolution */ 
  PointingResolution pr; 
  p.point->computePointingResolution(sum,pol, peak, &pr);  
  contribution.clear(); 
  // std::cout <<"\tdphi="<< pr.getdPhi()<<" dTheta="<< pr.getdTheta()<<std::endl;


  if ( pr.getdPhi() > p.max_dphi || pr.getdTheta() > p.max_dtheta) 
  {
    if (p.verbosity > 0) 
    {
      printf("Rejecting %d:%d:%d due to failing pointing resolution cut (dphi: %g, dtheta: %g\n", sum->eventNumber, (int) pol, peak, pr.getdPhi(), pr.getdTheta()); 
    }

    return 0 ;
  }

  double phi0 = sum->peak[pol][peak].phi; 

  double inv_two_pi_sqrt_det = get_inv_two_pi_sqrt_det(pr.getdPhi(), pr.getdTheta(), pr.getCorr()); 
  double min_p = dist2dens(maxDistance(), inv_two_pi_sqrt_det); 

  std::vector<int> used ( segmentationScheme()->NSegments()); //each seg has an integer "used".
  UsefulAdu5Pat pat(gps); 

  if (p.projection == Params::BACKWARD) 
  {
    //start with guess
    const AnitaEventSummary::PointingHypothesis *pk = &sum->peak[pol][peak];  
    PayloadParameters guess;  
    int status =  PayloadParameters::findSourceOnContinent(pk->theta,pk->phi,gps, &guess, p.refract, p.collision_detection ? p.collision_params.dx : 0); 
    if (p.verbosity > 2) 
    {
      guess.source.to(AntarcticCoord::STEREOGRAPHIC); 
      printf("status =%d , loc = %g %g %g\n", status, guess.source.x, guess.source.y, guess.source.z); 
    }
    /* set up vector of segments to check */
    size_t nchecked = 0; 
    std::vector<int> segs_to_check; 
    segs_to_check.reserve(100);  // a plausible number 100 segments to check
    int guess_seg = segmentationScheme()->getSegmentIndex(guess.source); //guess the initial seg
    if (guess_seg < 0) return inv_two_pi_sqrt_det; 
    used[guess_seg] = 1; 
    segs_to_check.push_back(guess_seg); // put the initial guess seg into segments to check.

    AntarcticCoord pos0 = segmentationScheme()->getSegmentCenter(segs_to_check[0]); //position is the first segment
    pos0.to(AntarcticCoord::WGS84); 

    if (p.verbosity > 2) 
    {
      printf("center position: %g %g %g\n", pos0.x, pos0.y, pos0.z); 
    }
    int nsamples = p.backwards_params.num_samples_per_bin; 
    /* set up vectors for sample phis /thetas/ densities */ 

    //Figure out the enhancement steps; 
    std::vector<int> enhancement_steps (1+p.backwards_params.max_enhance); 
    enhancement_steps[0] = nsamples; 
    for (int e = 1; e <= p.backwards_params.max_enhance; e++)
    {
        enhancement_steps[e] = pow( ceil (sqrt(enhancement_steps[e-1] * 2)),2); 
    }

    int max_samples = enhancement_steps[p.backwards_params.max_enhance]; 




    std::vector<AntarcticCoord> samples(max_samples); 
    std::vector<double> dens(max_samples);; 
    std::vector<double> phis(max_samples);; 
    std::vector<double> thetas(max_samples);; 
    std::vector<double> xs(max_samples);; 
    std::vector<double> ys(max_samples);; 
    std::vector<bool> occluded(max_samples); 
    AntarcticCoord pos ; //position of the cur segment


    /** Loop over segments we need to check to see if p > cutoff */
    //initially, there is only the first segment to check. But in this while loop, more neighbour will be added into this segments to check list. 
    //So more and more segs will be checked until the no segsments has probability larger than a threshold.
    while (nchecked < segs_to_check.size())
    {
      int seg = segs_to_check[nchecked++]; 


      bool done_with_this_segment = false; 
      //if we are at a steep angle, we should enhance anyway 
      int enhance_factor =  TMath::Min(p.backwards_params.max_enhance ,  int(sum->peak[pol][peak].theta / 10)); 
      int noccluded = 0; 
      double seg_p = 0;
      double max_dens = 0; 

      while(!done_with_this_segment) 
      {
        //ENHANCE 
        // we want to find the smallest perfect square at least twice as many nsamples
        if (enhance_factor) 
        {
          nsamples = enhancement_steps[enhance_factor]; 
          if (p.verbosity > 0) printf("ENHANCE! To %d\n", nsamples); 
        }
      
        // segment into a bunch of positions
        segmentationScheme()->sampleSegment(seg, nsamples, &samples[0], p.backwards_params.random_samples);  // do we want to randomize? i dunno. I'll decide later. 
        
        // loop over the samples, 
        //  we want to check:
        //     - can this sample probably see ANITA? 
        //     - what are the coordinates of this sample in ANITA's frame? 
        //
        for (int i = 0; i < nsamples; i++)
        {
          //This computes the payload in source coords and vice versa
          PayloadParameters pp(gps, samples[i], p.refract); 


          pos = samples[i].as(AntarcticCoord::STEREOGRAPHIC); 

          if (p.verbosity > 3) 
          {
            printf("   sample %d position: %g %g %g\n", i, pos.x, pos.y, pos.z); 
            printf("      delta phi: %g\n",pp.source_phi - sum->peak[pol][peak].phi); 
            printf("      delta el: %g\n",pp.source_theta - sum->peak[pol][peak].theta); 
            printf("      payload el: %g %g\n", pp.payload_el, p.backwards_params.el_cutoff); 
          }
          //payload is below horizon.  or collides . So anita woundnt see the sample. So set occluded to 1.
          AntarcticCoord collid, collid_exit; 
          if (pp.payload_el < p.backwards_params.el_cutoff || ( p.collision_detection && pp.checkForCollision(p.collision_params.dx,&collid,&collid_exit, p.dataset, p.collision_params.grace))) 
          {
            if (pp.payload_el >=p.backwards_params.el_cutoff)  //colides
            {
               collid.to(AntarcticCoord::STEREOGRAPHIC); 
               if (p.verbosity > 3) 
               {
                 printf("Collision at (%g %g %g), surface at %g\n", collid.x,collid.y,collid.z, RampdemReader::SurfaceAboveGeoidEN(collid.x,collid.y, p.dataset)); 
               }


               // we want to make sure that whatever segment is occluding is is considered, if it isn't already. So let's project to continent from payload and ensure we have that segment already 
               int potential_seg = p.seg->getSegmentIndex(collid_exit); 
               if (!used[potential_seg]) 
               {
                   used[potential_seg] = segs_to_check.size(); 
                   segs_to_check.push_back(potential_seg); 
               }
              
            }
            occluded[i] = true; 
          }
          else occluded[i] = false; 

          /* add to phis and thetas for this segment */

          phis[i]=pp.source_phi; 

          // Have to do something funny with phi if it's too different from initial phi 
          while(phis[i] - phi0 > 180) phis[i] -=360; 
          while(phis[i] - phi0 < -180) phis[i] +=360; 

          thetas[i]=pp.source_theta; 
          xs[i] = pos.x; 
          ys[i] = pos.y; 
        }
      

        /* compute the probabilities for each set of angles */ 
        //Get the Prob density for each sample in this segment.
        pr.computeProbabilityDensity(nsamples, &phis[0], &thetas[0], &dens[0]); 

        /** Set occluded's sample's prob density to 0
         *
         * Why not just get rid of them entirely you ask? I think the zero's might be important
         * for calculating the integral properly. But maybe they're not. I'll look into it. 
         *
         * */ 
        noccluded = 0; 
        max_dens = 0;
        for (int i =0; i < nsamples; i++)
        {
          if (occluded[i])
          {
            dens[i] = 0; 
            noccluded++; 
          }
          else if (dens[i] > max_dens)
          {
            //get the maximum density for this segment.
            max_dens = dens[i]; 
          }
        }



        /* Now we want to compute the integral  */

        /** write out triangles for debugging if wanted*/ 
        if (debugfile)
        {
          if (!debugfile->Get("triangles_payload")) debugfile->mkdir("triangles_payload"); 
          debugfile->cd("triangles_payload"); 
          TGraph2D g2d (nsamples,&phis[0], &thetas[0], &dens[0]); 
          g2d.Write(TString::Format("g%d",seg)); 

          if (!debugfile->Get("triangles_continent")) debugfile->mkdir("triangles_continent"); 
          debugfile->cd("triangles_continent"); 

          TGraph2D g2d_continent (nsamples,&xs[0], &ys[0], &dens[0]); 
          g2d_continent.Write(TString::Format("g%d",seg)); 

        }

        bool enhance_flag = false;
        double scale = 0; 
        #ifdef HAVE_DELAUNAY
        // using the n samples to integrate the probabilty for this segment. The math tool is delanuay2d.
        ROOT::Math::Delaunay2D del (nsamples,&phis[0], &thetas[0], &dens[0]); 
        del.FindAllTriangles(); 
        // the delaunay triangulation stupidly normalizes... have to unnormalize it
        double max_phi = -360000; 
        double min_phi = 360000; 
        double max_theta = -360000; 
        double min_theta = 360000; 
        for (int i = 0; i < nsamples; i++) 
        {
          if (phis[i] > max_phi) max_phi = phis[i]; 
          if (thetas[i] > max_theta) max_theta = thetas[i]; 
          if (phis[i] < min_phi) min_phi = phis[i]; 
          if (thetas[i] < min_theta) min_theta = thetas[i]; 
        }

        scale=  (max_phi-min_phi) * (max_theta-min_theta)  / ((del.XMax()-del.XMin())*(del.YMax()-del.YMin())); 

        double sum = 0; 
        for (std::vector<ROOT::Math::Delaunay2D::Triangle>::const_iterator it = del.begin(); it!=del.end(); it++)
        {
          const ROOT::Math::Delaunay2D::Triangle & tri = *it; 

          double area = 0.5 *scale* fabs(( tri.x[0] - tri.x[2]) * (tri.y[1] - tri.y[0]) - (tri.x[0] - tri.x[1]) * (tri.y[2] - tri.y[0])); 
          for (int j = 0; j < 3; j++) 
          {
            double darea =area /3. *  dens[tri.idx[j]];
            if (p.backwards_params.enhance_threshold < darea && enhance_factor < p.backwards_params.max_enhance) 
            {
              enhance_flag = true; 
              break; 
            }
            sum += darea; 
          }
        }
        #else
        double sum =-1; 
        fprintf(stderr,"ROOT 5 not currently supported in Probability Map due to lack of ROOT/Math/Delaunay2D.h. Will be fixed someday if necessary. \n"); 
        #endif
        //a segment is an area, intergrate over this area can give us the probabilty of signal coming from this segmentl.
        seg_p = sum; 
        // printf("%d %g %d\n",seg, seg_p, noccluded); 
 

        if (p.verbosity > 2) 
        {

          printf("p(seg) %d: = %g\n",  seg, seg_p); 

        }
        //seg_p is the intergral 
        if (seg_p > 1) 
        {
          // will stay in the while loop if the enhancement is not reached.
          if ( enhance_factor < p.backwards_params.max_enhance) 
          {
            enhance_flag = true; 
          }
          else
          {
            printf("OOPS segment(%d) scale(%g) seg_p=%g>1, just set it to 1.\n", seg, scale, seg_p); 
            seg_p = 1; 
          }
        }

        done_with_this_segment = !enhance_flag; 
        enhance_factor++; //TODO, we can enhance more than once depending on how much above threshold we are 
      }
      // is an vertor store the segment number and fraction of sample that are occluded.
      if (occlusion) 
      {
          occlusion->push_back(std::pair<int,double> ( seg, double(noccluded)/nsamples)); 
      }


      pos = segmentationScheme()->getSegmentCenter(seg);
      pos0.to(AntarcticCoord::STEREOGRAPHIC); 
      pos.to(AntarcticCoord::STEREOGRAPHIC); 
      
      /* if max density is above min_p, add the neighbors of this segment */ 
      // This is done recursively so the neighbour segments(which has max_dens> min_p) will be added into the contribution vector.
      //min_p is calculated from the max distance, like 20 sigma. So it is the upper limit for a point.
      // std::cout << pos.x<< " " << pos0.x << " "<< pos.y << " "<< pos0.y << " " << sqrt((pos.x - pos0.x)*(pos.x - pos0.x) + (pos.y - pos0.y)*(pos.y - pos0.y))<<  std::endl;
      if (max_dens >= min_p or sqrt((pos.x - pos0.x)*(pos.x - pos0.x) + (pos.y - pos0.y)*(pos.y - pos0.y)) <= p.radius)
      {
        //moved these two line inside the if condition. Because we only record the segment if it's max density > min_p or  distance from base to projection position < radius,
        //is an vertor store the current segment and its prob density integral over the area.
        contribution.push_back(std::pair<int,double>(seg,seg_p));
        //is an vector store the current segment number and it maximun den sity 
        if(max_densities) max_densities->push_back(std::pair<int,double>(seg, max_dens)); 

        std::vector<int> new_neighbors;
        segmentationScheme()->getNeighbors(seg, &new_neighbors); // find the 8 neighbor seg around this seg. The return vector is in new_neighbors.
        for (size_t j = 0; j < new_neighbors.size(); j++)
        {
          int new_seg = new_neighbors[j];
          //each time a new neighbor segment need to be check, it should be also note as used when it is not used.. 
          //Used also works as the count of the segments to check list. So it will increase as segstocheck size increase.
          if (!used[new_seg])
          {
            used[new_seg] = segs_to_check.size(); 
            segs_to_check.push_back(new_seg); 
          }
        }
      }
    }

    //finally, loop over any bases (if we were asked to do so) 
    if (base_contribution) 
    {
      // std::cout<< min_p << std::endl;
      //now loop over all the bases 
      int nbases = BaseList::getNumAbstractBases(); 
      for (int ibase = 0; ibase < nbases; ibase++)
      {

        const BaseList::abstract_base & base = BaseList::getAbstractBase(ibase); 
        AntarcticCoord base_pos= base.getPosition(gps->realTime);

        //TODO: can veto most bases early probalby 
        //give the gps and base positon, easy to figure out all the geom between payload and base.
        PayloadParameters pp (gps,base_pos, p.refract);
        // if this base can not see the payload over horizon , will continue to look for next base. 
        // if (pp.payload_el < p.backwards_params.el_cutoff || ( p.collision_detection && pp.checkForCollision(p.collision_params.dx,0,0, p.dataset, p.collision_params.grace))) continue ; 
        if (pp.payload_el < p.backwards_params.el_cutoff) {
          continue;
        }
        // when a base is in view, compute its prob density at this base point.
        double base_phi = pp.source_phi ; 
        if (base_phi - phi0 > 180) base_phi-=360; 
        if (base_phi - phi0 < -180) base_phi+=360; 
        double dens = pr.computeProbabilityDensity( base_phi, pp.source_theta); 
        // std::cout<< ibase << " " << dens << " "  << std::endl;
        //this is vector that record the base id and its prob density.
        base_pos.to(AntarcticCoord::STEREOGRAPHIC); 

        if(dens> min_p or sqrt((base_pos.x - pos0.x)*(base_pos.x - pos0.x) + (base_pos.y - pos0.y)*(base_pos.y - pos0.y)) <= 2 * p.radius){
          // only when dens larger than the min_p or distance from base to projection position < 2 * radius, it will record the base.
          // the size of base_contribution will tell us the number of bases that this event is overlapping with.
          //overlapping is defined the same as the p.maximum_distance(ie, how many sigma)
          base_contribution->push_back(std::pair<int,double> (ibase, dens));
        }
      }

    }

  }

  // this MC method dose work anymore, need to remove.
  else if (p.projection == Params::MC) 
  {

    std::vector<long long> counts(segmentationScheme()->NSegments()); 
    for (long long i = 0; i < p.mc_params.n; i++) 
    {

      double phi, theta; 
      pr.random(&phi,&theta); 
      //      double lon,lat,alt; 
      //        printf("%lld %f %f\n",i, phi,theta); 
      //
      //      int success = pat.getSourceLonAndLatAtAlt(phi * TMath::DegToRad(), theta * TMath::DegToRad(), lon,lat,alt); 

      PayloadParameters pp; 

      int status = PayloadParameters::findSourceOnContinent(theta,phi,gps, &pp, p.refract, p.collision_detection ? p.collision_params.dx :0); 
      if (status == 1) 
      {
      //AntarcticCoord c(AntarcticCoord::WGS84,lat,lon,alt); 
         int seg = p.seg->getSegmentIndex(pp.source); 
         counts[seg]++; 
      //printf("%d %lld\n", seg, counts[seg]); 
      }
    }

    for (size_t i = 0; i < counts.size(); i++)
    {
      if (counts[i])
      {
        contribution.push_back(std::pair<int,double>(i, double(counts[i]) / p.mc_params.n)); 
        if(occlusion) occlusion->push_back(std::pair<int,double>(i,0)); 
      }
    }

  }


  
  
  return inv_two_pi_sqrt_det; 
}

// int UCorrelator::ProbabilityMap::groupAdjacent(const double * vals_to_group, std::vector<std::vector<int> > * groups, double * counts, std::vector<double>  * dist, double val_threshold) const
// {
//   //vals_to_group are map from seg index to seg's ps. All the wgt from one event should clustered and sum to 1. 
//   std::vector<bool> consumed (p.seg->NSegments()); 
//   int ngroups = 0; 


//   //clean up input 
//   if (groups) groups->clear(); 
//   if (dist) dist->clear(); 
//   if (counts) memset(counts, 0, sizeof(double) * p.seg->NSegments()); 

//   for (int i = 0; i < p.seg->NSegments(); i++) 
//   {
//     //only clustering segments when the segment's weight larger than a threshold.
//     if (vals_to_group[i]<=val_threshold){
//       consumed[i] = true;
//       continue;
//     } 
//      if (consumed[i]) continue; 

//      ngroups++; 
//      //group is a vector of int, put inital seg index in it at first, then add the neighbor seg index recursively.
//      std::vector<int> group;
//      group.push_back(i); 
//      consumed[i] = true; 
//      int n_in_group = 0; 
//      //i start with 0. group_sum initially will have the first segment's ps value.
//      double group_sum = vals_to_group[i]; 
//      while (n_in_group < (int)  group.size()) 
//      {
//        int seg = group[n_in_group++]; 

//        //see if we need to add any neighbors 

//        std::vector<int> neighbors; 
//        segmentationScheme()->getNeighbors(seg, &neighbors); 
//        for (size_t i = 0; i < neighbors.size(); i++) 
//        {
//         // only cluster the neighbour segment when their ps larger than threshold and not consumed before.
//          if (vals_to_group[neighbors[i]]>val_threshold && !consumed[neighbors[i]])
//          {
//            consumed[neighbors[i]] = true;
//            group.push_back(neighbors[i]);
//            //group_sum added the neighbour segs' ps togother
//            group_sum += vals_to_group[neighbors[i]]; 
//          }
//        }
//      }
//      //group_sum is the number of events(prob sums together = number of events) in this cluster. It is calculated by clustering the segments and summing the weights. 
//      //dist is a vector of length = how many clusters.
//      if (dist) dist->push_back(group_sum); 

//      if (counts) 
//      {
//       // loop through the current cluster's segment index i.
//        for (size_t i = 0; i < group.size(); i++) 
//        {
//         //group[i] is one segment's index number
//         //counts is a map from segment index to a cluster's prob sums(the cluster that this segment is in.)
//          counts[group[i]] = group_sum; 
//        }
//      }

//      if (groups) groups->push_back(group); 
//   }
  
//   return ngroups; 

// }
int UCorrelator::ProbabilityMap::countBasesInThisSegment(int seg) const{
  int countBases = 0;
  //now loop over all the bases 
  int nbases = BaseList::getNumBases(); 
  for (int ibase = 0; ibase < nbases; ibase++)
  {

    const BaseList::abstract_base & base = BaseList::getBase(ibase);
    // time is 0 and irrelavent to stationary bases. 
    AntarcticCoord base_pos= base.getPosition(0);
    if (segmentationScheme()->getSegmentIndex(base_pos) == seg){
      countBases += 1;
    }
  }
  return countBases;
}
int UCorrelator::ProbabilityMap::doClustering(double threshold)
{
  //ps are map from seg index to seg's ps. All the wgt from one event should clustered and sum to 1. 
  std::vector<bool> consumed (p.seg->NSegments()); 
  int nClusters = 0; 

  for (int i = 0; i < p.seg->NSegments(); i++) 
  {
    //only clustering segments when the segment's weight larger than a threshold.
    if (uniform_ps[i]<=threshold){
      consumed[i] = true;
      continue;
    } 
     if (consumed[i]) continue; 
     

     nClusters++; 
     //Cluster is a vector of int, put inital seg index in it at first, then add the neighbor seg index recursively.
     std::vector<int> Cluster;
     Cluster.push_back(i); 
     consumed[i] = true; 
     int n_in_Cluster = 0; 
     //i start with 0. cluster_sum initially will have the first segment's uniform_ps value.
     double cluster_sum = uniform_ps[i];
     double cluster_sum_with_base = uniform_ps_with_base[i];
     double cluster_sum_without_base = uniform_ps_without_base[i];

     while (n_in_Cluster < (int)  Cluster.size()) 
     {
       int seg = Cluster[n_in_Cluster++]; 

       //see if we need to add any neighbors 

       std::vector<int> neighbors; 
       segmentationScheme()->getNeighbors(seg, &neighbors); 
       for (size_t j = 0; j < neighbors.size(); j++) 
       {
        // only cluster the neighbour segment when their uniform_ps larger than threshold and not consumed before.
         if (uniform_ps[neighbors[j]]>threshold && !consumed[neighbors[j]])
         {
           consumed[neighbors[j]] = true;
           Cluster.push_back(neighbors[j]);
           //cluster_sum added the neighbour segs' uniform_ps togother
           cluster_sum += uniform_ps[neighbors[j]];
           cluster_sum_with_base += uniform_ps_with_base[neighbors[j]];
           cluster_sum_without_base += uniform_ps_without_base[neighbors[j]];
         }
       }
     }
     //finish one cluster. 

     //cluster_sum is the number of events(prob sums together = number of events) in this cluster. It is calculated by clustering the segments and summing the weights. 
     //clusterSizes is a vector of length = how many clusters.
     clusterSizes.push_back(cluster_sum); 
     clusterSizes_with_base.push_back(cluster_sum_with_base); 
     clusterSizes_without_base.push_back(cluster_sum_without_base); 

    // loop through the current cluster's segment index i.
     for (size_t i = 0; i < Cluster.size(); i++) 
     {
      //Cluster[i] is one segment's index number
      //mapOfClusterSizes is a map from segment index to a cluster's prob sums(the cluster that this segment is in.)
       mapOfClusterSizes[Cluster[i]] = cluster_sum; 
       mapOfClusterIndexs[Cluster[i]] = nClusters; 
     }
  }
  
  return nClusters; 

}

void  UCorrelator::ProbabilityMap::evaluateEvent(double & indexOfCluster, double & sizeOfCluster, double theta,double phi,const Adu5Pat * gps){
  PayloadParameters guess;  
  int status =  PayloadParameters::findSourceOnContinent(theta,phi,gps, &guess, p.refract, p.collision_detection ? p.collision_params.dx : 0);
  // if (status == 0){
  //   // over horizon
  //   return -1;
  // }else{
    guess.source.to(AntarcticCoord::STEREOGRAPHIC); 
    int seg = segmentationScheme()->getSegmentIndex(guess.source);
    indexOfCluster = mapOfClusterIndexs[seg];
    sizeOfCluster = mapOfClusterSizes[seg];
    return;
  // } 
}

std::pair<int, int> UCorrelator::ProbabilityMap::showClusters(int draw, bool blind, const char * option) const
{
  //hard coded to blind, comment out to unblind
  blind = false;

  const int maxes[] = {1,2,3,4,5,6,7,8,9,10,20,50,100,(int) 100e6}; 
  const int mins[] =  {1,2,3,4,5,6,7,8,9,10,11,21,51,101}; 
  int n_rows = sizeof(maxes) / sizeof(*maxes); 

  int n_clusters[n_rows]; 
  int n_clusters_near_base[n_rows]; 
  int n_clusters_not_base[n_rows]; 
  float n_clusters_weighted[n_rows]; 
  memset(n_clusters,0, n_rows * sizeof(int)); 
  memset(n_clusters_near_base,0, n_rows * sizeof(int)); 
  memset(n_clusters_not_base,0, n_rows * sizeof(int)); 
  memset(n_clusters_weighted,0, n_rows * sizeof(float)); 
  for (size_t i = 0; i < clusterSizes.size(); i++) 
  {
    int n = round(clusterSizes[i]); 
    int n_nearBase =round(clusterSizes_with_base[i]);
    int n_notnearBase = round(clusterSizes_without_base[i]); 
    if(blind and n==1 and n_notnearBase==1){
      // if the case is what we should be blind, then skip the output.
      ;
    }else{
      std::cout <<n << " \t"<< n_nearBase<< " \t" <<n_notnearBase <<  std::endl;
    }

    float fractionOfEventsNearBase = n_nearBase/float(n);

    for (int row = 0; row < n_rows; row++)
    {
       if (n <= maxes[row] && n >= mins[row])
       {
         n_clusters[row]++;
         if(fractionOfEventsNearBase!=0){
          n_clusters_near_base[row]++;
         }else{
          n_clusters_not_base[row]++;
         }
         n_clusters_weighted[row] +=  fractionOfEventsNearBase;
         break; 
       }
    }
  }
  // if we need to be blind to Signal box, just set this element to 0.
  if(blind){
    n_clusters[0] -= n_clusters_not_base[0];
    n_clusters_not_base[0] = 0;
  }

  if (draw){
    printf("================================================================================================================\n"); 
    printf("|Cluster Size\t|N of Clusters\t|N of weighted Clusters\t|N of Clusters nearBase\t|N of Clusters not near Base\t|\n");
    printf("----------------------------------------------------------------------------------------------------------------\n");
    for (int row = 0; row < n_rows; row++) 
    {
      if(mins[row] == maxes[row]){
        printf ("|     %d      ", mins[row]);
      }
      else if(row < n_rows-1)
        printf ("|  %d - %d   ", mins[row], maxes[row]); 
      else
        printf("|  %d +       ",mins[row]);
      printf("\t|      %d     \t|     %f    \t|           %d          \t|              %d           \t|\n", 
                      n_clusters[row], n_clusters_weighted[row], n_clusters_near_base[row],n_clusters_not_base[row]);

       
    }
    printf("----------------------------------------------------------------------------------------------------------------\n");  
  }

  
  std::vector <double> temp_map(p.seg->NSegments(),0);

  if(draw == 1){
    temp_map = ps_norm;
  }else if(draw == 2){
    temp_map = mapOfClusterSizes;
  }else if(draw == 3){
    temp_map = mapOfClusterIndexs;
  }else if(draw == 4){
    temp_map = uniform_ps_weighted_by_base;
  }else if(draw == 5){
    temp_map = uniform_ps_with_base;
  }else if(draw == 6){
    temp_map = uniform_ps_without_base;
  }
  // blinded to the map unknown singlets
  if(blind){
    for (int i =0; i< p.seg->NSegments(); i++){
      if(round(mapOfClusterSizes[i]) == 1 and uniform_ps_without_base[i] != 0){
        temp_map[i] = 0;
      }
    }
  }
  segmentationScheme()->Draw(option,&temp_map[0]);


  // calculate the total number of clusters near base
  double sumNumOfClusterNearBase = 0;
  for (int row =0; row < n_rows; row++){
    sumNumOfClusterNearBase+= n_clusters_near_base[row];
  }



  // TCanvas * canvas = new TCanvas;
  // canvas->Divide(1,1); 
  // canvas->cd(1);
  // TH1F *h1 = new TH1F("h1", "h1 title", 20, 0, 5);
  double sumFractionOfEventsNearBase = 0;
  int nSinglets_knownBase = 0;
  int nSinglets_unkownBase = 0;
  for (int i = 0; i < clusterSizes.size(); i++) 
  {
    int n = round(clusterSizes[i]); 
    int n_nearBase =round(clusterSizes_with_base[i]);
    int n_notnearBase = round(clusterSizes_without_base[i]); 
    float fractionOfEventsNearBase = n_nearBase/float(n);
    sumFractionOfEventsNearBase +=fractionOfEventsNearBase;
    // std::cout<< "\tcluster i="<< i << "\t #events = "<< n <<  "\t #eventNearBase = "<< n_nearBase<< "\t fractionOfEventsNearBase = "<< fractionOfEventsNearBase << std::endl;
    if(n == 1 and n_nearBase == 1){
      nSinglets_knownBase++;
    }else if(n == 1 and n_notnearBase == 1){
      nSinglets_unkownBase++;
    }
  }
  if (blind){
    nSinglets_unkownBase = 0;
  }

  // std::cout<< countNofClusterWithBase<< " \t" << countNofClusterWithoutBase<< " \t"<< countNofEventsWithBase<< " \t"<< countNofEventsWithoutBase<< " \t"<< countNofBase << " \t"<< countNofUnkownBase<< std::endl;

  //background estimate:
  int A=0, B =0;
  int C = n_clusters_not_base[0];
  int D = n_clusters_near_base[0];
  float B1 = 0;
  for (int row = 1; row < 5; row++){
    A += n_clusters_not_base[row];
    B += n_clusters_near_base[row];
    B1 += n_clusters_weighted[row];
  }
  float C0 = float(D)*float(A)/float(B);
  float C1 = float(D)*float(A)/float(B1);
  std::cout<< "Anthropogenic background(ABCD mehtod)="<< C0<< "\t Anthropogenic background(Weighted ABCD mehtod)="<< C1<< std::endl;
  std::cout<<"\n\n" <<std::endl;
  std::cout<<sumNumOfClusterNearBase<< " \t"<<sumFractionOfEventsNearBase<< " \t"<< A<< " \t"<< B << " \t"<< C<<" \t"<< D<< " \t"<< C0 << " \t"<< C1 << std::endl;
  // TGraph *gr1 = new TGraph(N, x, y);
  // h1->Draw("colz");
  return std::make_pair(nSinglets_knownBase, nSinglets_unkownBase);
  // return n_clusters_not_base[0]; 

}


double UCorrelator::ProbabilityMap::getProbSumsIntegral(bool norm) const
{
  const double * V = getProbSums(norm); 
  int N = segmentationScheme()->NSegments(); 
  double sum = 0; 
  for (int i = 0; i < N; i++) 
  {
    sum += V[i]; 
  }

  return sum; 
}



