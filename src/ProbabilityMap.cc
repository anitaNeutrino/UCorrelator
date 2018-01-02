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
    ps_without_base(NLevels(), std::vector<double> (p.seg->NSegments(),0)), 
    ps_norm(p.seg->NSegments(),0), 
    ps_norm_without_base(NLevels(), std::vector<double> (p.seg->NSegments(),0)), 
    max1_ps(p.seg->NSegments(),0), 
    max1_ps_norm(p.seg->NSegments(),0), 
    max2_ps(p.seg->NSegments(),0), 
    max2_ps_norm(p.seg->NSegments(),0), 
    sqrt_ps(p.seg->NSegments(),0), 
    sqrt_ps_without_base(NLevels(), std::vector<double> (p.seg->NSegments(),0)), 
    sqrt_ps_norm(p.seg->NSegments(),0), 
    sqrt_ps_norm_without_base(NLevels(), std::vector<double> (p.seg->NSegments(),0)), 
    fraction_occluded(p.seg->NSegments(), 0), 
    n_above_level(NLevels(), std::vector<int>(p.seg->NSegments(),0)),
    n_above_level_norm(NLevels(), std::vector<int>(p.seg->NSegments(),0)),
    wgt_above_level(NLevels(), std::vector<double>(p.seg->NSegments(),0)),
    wgt_above_level_norm(NLevels(), std::vector<double>(p.seg->NSegments(),0)),
    n_above_level_without_base(NLevels(), std::vector<int>(p.seg->NSegments(),0)),
    n_above_level_without_base_norm(NLevels(), std::vector<int>(p.seg->NSegments(),0)),
    wgt_above_level_without_base(NLevels(), std::vector<double>(p.seg->NSegments(),0)),
    wgt_above_level_without_base_norm(NLevels(), std::vector<double>(p.seg->NSegments(),0)),
    base_n_above_level(NLevels(), std::vector<int>(BaseList::getNumBases() + BaseList::getNumPaths(), 0)), 
    base_n_above_level_norm(NLevels(), std::vector<int>(BaseList::getNumBases() + BaseList::getNumPaths(), 0)), 
    base_sums(BaseList::getNumBases() + BaseList::getNumPaths()) , 
    base_sums_norm(BaseList::getNumBases() + BaseList::getNumPaths()) 
{


}


int UCorrelator::ProbabilityMap::add(const AnitaEventSummary * sum, const Adu5Pat * pat, AnitaPol::AnitaPol_t pol,
                                     int peak, double weight, TFile * debugfile) 
{

  std::vector<std::pair<int,double> > segments_to_fill; 
  std::vector<std::pair<int,double> > base_ps_to_fill; 
  std::vector<std::pair<int,double> > occluded_to_fill; 
  std::vector<std::pair<int,double> > max_densities; ; 

  double inv_two_pi_sqrt_det = computeContributions(sum,pat,pol,peak,segments_to_fill, &base_ps_to_fill, &occluded_to_fill, &max_densities, debugfile); 

  if (!inv_two_pi_sqrt_det) return 0; 

  std::vector<double> levels_p(NLevels()); 

  for (int i = 0; i < int(NLevels()); i++)
  {
    levels_p[i] = dist2dens(p.level_thresholds[i], inv_two_pi_sqrt_det); 
    if (p.verbosity > 2) 
    {
      printf("mahalanobis distance: %g, probability density: %g\n", p.level_thresholds[i], levels_p[i]); 
    }
  }
  
  TLockGuard lock(&m); 

  int incr = weight > 0 ? 1 : weight < 0 ? -1 : 0; 
  int Nbases = base_ps_to_fill.size(); 
  int min_base_level = NLevels(); 
  int min_base_level_norm = NLevels(); 

  int NFilled = segments_to_fill.size(); 
  double norm = 0; 
  for (int i = 0; i < NFilled; i++)
  {
    norm += segments_to_fill[i].second; 
  }

  double invnorm = norm < p.min_p_on_continent ? 0 : 1./norm;
  if (p.verbosity > 2) printf("invnorm: %g\n", invnorm); 

  for (int i = 0; i < Nbases; i++)
  {
    int ibase = base_ps_to_fill[i].first; 
    double dens_base = base_ps_to_fill[i].second * weight; 
    base_sums[ibase] += dens_base;
    base_sums_norm[ibase] += dens_base * invnorm; 

    for (int j = 0; j < (int) NLevels(); j++)
    {
      if (dens_base > levels_p[j])
      {
        base_n_above_level[j][ibase]+=incr; 
        if (min_base_level > j) min_base_level = j; 
      }

      if (dens_base * invnorm > levels_p[j])
      {
        base_n_above_level_norm[j][ibase]+=incr; 
        if (min_base_level_norm > j) min_base_level_norm = j; 
      }
    }
  }




  std::vector<int> this_NAboveLevel(NLevels()); 
  std::vector<int> this_NAboveLevel_norm(NLevels()); 

  for (int i = 0; i < NFilled; i++)
  {
    double max_dens_seg = max_densities[i].second * weight; 
    double max_dens_seg_norm = max_densities[i].second * invnorm * weight; 
    for (int j = 0; j < (int) NLevels(); j++) 
    {
      if (max_dens_seg > levels_p[j])
      {
        this_NAboveLevel[j]+=incr; 
      }

      if (max_dens_seg_norm > levels_p[j])
      {
        this_NAboveLevel_norm[j]+=incr; 
      }
    }
  }


  for (int i = 0; i < NFilled; i++)
  {
    int iseg = segments_to_fill[i].first; 
    double pseg = segments_to_fill[i].second * weight; 
    double max_dens_seg = max_densities[i].second * weight; 
    double max_dens_seg_norm = max_dens_seg * invnorm; 

    double pseg_norm = pseg * invnorm; 

    if (pseg > max1_ps[iseg] ) 
    {
      max2_ps[iseg] = max1_ps[iseg];
      max1_ps[iseg ] = pseg ; 
    }
    else if (pseg > max2_ps[iseg]) 
    {
      max2_ps[iseg] = pseg ;
    }

    if (pseg_norm > max1_ps_norm[iseg] ) 
    {
      max2_ps_norm[iseg] = max1_ps_norm[iseg];
      max1_ps_norm[iseg ] = pseg_norm ; 
    }
    else if (pseg  > max2_ps_norm[iseg]) 
    {
      max2_ps_norm[iseg] = pseg_norm ;
    }


    ps[iseg] += pseg; 
    ps_norm[iseg] += (pseg_norm); 
    sqrt_ps[iseg]  += weight < 0 ? -sqrt(-pseg) : sqrt(pseg); 
    sqrt_ps_norm[iseg] +=  weight < 0 ? -sqrt(-pseg_norm) : sqrt(pseg_norm); 

    for (int level = 0; level < min_base_level; level++) 
    {
      ps_without_base[level][iseg] += (pseg);
      sqrt_ps_without_base[level][iseg] += weight < 0 ? -sqrt(-pseg) : sqrt(pseg); 
    }
    for (int level = 0; level < min_base_level_norm; level++)
    {
      ps_norm_without_base[level][iseg] += (pseg_norm);
      sqrt_ps_norm_without_base[level][iseg] += weight < 0 ? -sqrt(-pseg_norm): sqrt(pseg_norm);
    }

    fraction_occluded[iseg] += occluded_to_fill[i].second; 

    for (int j = 0; j < (int) NLevels(); j++)
    {
      if (max_dens_seg > levels_p[j]) 
      {
//        printf("Segment %d above level %d (%g) \n", iseg, j, pseg); 
        n_above_level[j][iseg]+=incr; 
        wgt_above_level[j][iseg]+= incr/ double(this_NAboveLevel[j]); 
        if (min_base_level > j)
        {
          n_above_level_without_base[j][iseg]+=incr; 
          wgt_above_level_without_base[j][iseg]+= incr/ double(this_NAboveLevel[j]); 
        }
      }
        
      if (max_dens_seg_norm > levels_p[j]) 
      {
//        printf("Segment %d above level %d (%g) \n", iseg, j, pseg); 
        n_above_level_norm[j][iseg]+=incr; 
        wgt_above_level_norm[j][iseg]+= incr/ double(this_NAboveLevel[j]); 
        if (min_base_level_norm > j)
        {
          n_above_level_without_base_norm[j][iseg]+=incr; 
          wgt_above_level_without_base_norm[j][iseg]+= incr/ double(this_NAboveLevel[j]); 
        }
      }
    }
  }


  return NFilled; 
}

int UCorrelator::ProbabilityMap::dumpNonZeroBases()  const 
{
  int n = 0;
  for (int i = 0; i < (int) getNBases(); i++)
  {
    if(base_sums[i] > 0) 
    {
      printf("%s (%d): %g\n", BaseList::getAbstractBase(i).getName(), i, base_sums[i]); 
      n++; 
    }
  }
 
  return n; 
}

double UCorrelator::ProbabilityMap::overlap(const AnitaEventSummary * sum , const Adu5Pat * pat,  AnitaPol::AnitaPol_t pol, int peak, 
                                            bool normalized, double weight,  std::vector<std::pair<int,double> > * bases,
                                            OverlapMode mode,  bool remove_self, std::vector<std::pair<int, double >> * segments, 
                                            std::vector<std::pair<int,double> > * max_dens, double * inv_two_pi_sqrt_det ) const
{
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
    if (!inv || norm < p.min_p_on_continent)
    {
      printf("below minimum: %g!\n", norm); 
      return -1; // there can be no overlap 
    }
    invnorm = 1./norm; 
  }

  if (inv_two_pi_sqrt_det) *inv_two_pi_sqrt_det = inv; 


  //pick the right thing to overlap with 
  const double * the_rest = mode == OVERLAP_SUM_SQRTS ? getProbSqrtSums(normalized) : 
                            mode == OVERLAP_SQRT_SUMS || mode == OVERLAP_SUMS ? getProbSums(normalized) : 
                            getProbMaxes(normalized); 

 for (int i =0; i< N; i++)
  {
    int seg = segs[i].first; 

    double p_this  =  weight*segs[i].second * invnorm; 

    double p_other = the_rest[seg]; 


    if (remove_self && mode == OVERLAP_SQRT_SUMS)
    {
      p_other -= p_this; 
    }

    if (remove_self && mode == OVERLAP_SUM_SQRTS) 
    {
      p_other -= sqrt(p_this);  //this might help with floating point precision to do first... maybe. 
    }

    if (remove_self && mode == OVERLAP_MAXES &&  p_this == p_other) p_other = normalized ? max2_ps_norm[seg] : max2_ps[seg]; 

    double danswer = (mode == OVERLAP_SUM_SQRTS ? sqrt(p_this) : p_this) * p_other; 
    
    if (mode == OVERLAP_SQRT_SUMS || mode == OVERLAP_MAXES) danswer = sqrt(danswer); 

    if (remove_self &&  mode == OVERLAP_SUMS)
    {
      danswer -= p_this; 
    }

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


  if (NLevels() != other.NLevels())
  {
    fprintf(stderr,"Cannot combine ProbabilityMap's with different numbers of levels: (ours %d, theirs: %d)\n", (int) NLevels(), (int) other.NLevels()); 
    return 1; 
  }

  int nwrong = 0; 
  
  for (int i = 0; i < (int) NLevels(); i++) 
  {
    if (getLevel(i) != other.getLevel(i))
    {
      fprintf(stderr,"  Mismatched level %d (ours=%g, theirs=%g)", i, getLevel(i), other.getLevel(i));
      nwrong++;
    }
  }

  if (nwrong > 0) 
  {
    fprintf(stderr,"%d Mismatched Levels. Cannot combine. Theoretically we could combine the levels that still make sense, but you'll have to implement that if you really want it.\n", nwrong); 
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
  

  //TODO we probably want to compare the pointing resolution models as well 

  //if we made it this far, all is good (except for what hasn't been implemented yet) !

  for (int i = 0; i < segmentationScheme()->NSegments(); i++) 
  {
    ps[i] += other.getProbSums()[i]; 
    ps_norm[i] += other.getProbSums(true)[i]; 
    sqrt_ps[i] += other.getProbSqrtSums()[i]; 
    sqrt_ps_norm[i] += other.getProbSqrtSums(true)[i]; 


    double other_max = TMath::Min(max1_ps[i], other.getProbMaxes()[i]); 
    max1_ps[i] = TMath::Max(max1_ps[i], other.getProbMaxes()[i]); 
    double other_max_norm = TMath::Min(max1_ps[i], other.getProbMaxes(true)[i]); 
    max1_ps_norm[i] = TMath::Max(max1_ps_norm[i], other.getProbMaxes(true)[i]); 

    max2_ps[i] = TMath::Max(other_max, TMath::Max(max2_ps[i], other.getProbSecondMaxes()[i])); 
    max2_ps_norm[i] = TMath::Max(other_max_norm, TMath::Max(max2_ps_norm[i], other.getProbSecondMaxes(true)[i])); 

    
    for (int j = 0; j < (int) NLevels(); j++)
    {
      ps_without_base[j][i] += other.getProbSumsWithoutBases(j)[i]; 
      ps_norm_without_base[j][i] += other.getProbSumsWithoutBases(j,true)[i]; 
      sqrt_ps_without_base[j][i] += other.getProbSqrtSumsWithoutBases(j)[i]; 
      sqrt_ps_norm_without_base[j][i] += other.getProbSqrtSumsWithoutBases(j,true)[i]; 
      n_above_level[j][i] += other.getNAboveLevel(j)[i]; 
      wgt_above_level[j][i] += other.getWgtAboveLevel(j)[i]; 
      n_above_level_without_base[j][i] += other.getNAboveLevelWithoutBases(j)[i]; 
      wgt_above_level_without_base[j][i] += other.getWgtAboveLevelWithoutBases(j)[i]; 
      n_above_level_without_base_norm[j][i] += other.getNAboveLevelWithoutBases(j,true)[i]; 
      wgt_above_level_without_base_norm[j][i] += other.getWgtAboveLevelWithoutBases(j,true)[i]; 
    }
  }

  if (combine_bases) 
  {
    for (size_t i = 0; i < getNBases(); i++) 
    {

      base_sums[i] += other.getBaseSums()[i]; 
      base_sums_norm[i] += other.getBaseSums(true)[i]; 
      for (int j = 0; j < (int) NLevels(); j++)
      {
        base_n_above_level[j][i] += other.getBaseNAboveLevel(j)[i]; 
        base_n_above_level_norm[j][i] += other.getBaseNAboveLevel(j,true)[i]; 
      }
    }
  }

  return 0; 
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

  std::vector<int> used ( segmentationScheme()->NSegments()); 
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
    segs_to_check.reserve(100);  // a plausible number
    int guess_seg = segmentationScheme()->getSegmentIndex(guess.source); 
    if (guess_seg < 0) return inv_two_pi_sqrt_det; 
    used[guess_seg] = 1; 
    segs_to_check.push_back(guess_seg); 

    AntarcticCoord pos = segmentationScheme()->getSegmentCenter(segs_to_check[0]); 
    pos.to(AntarcticCoord::WGS84); 

    if (p.verbosity > 2) 
    {
      printf("center position: %g %g %g\n", pos.x, pos.y, pos.z); 
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
    std::vector<bool> occluded(max_samples); 

    /** Loop over segments we need to check to see if p > cutoff */
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
//#pragma omp parallel for
        for (int i = 0; i < nsamples; i++)
        {
          //This computes the payload in source coords and vice versa
          PayloadParameters pp(gps, samples[i], p.refract); 


          pos = samples[i].as(AntarcticCoord::STEREOGRAPHIC); 
          if (p.verbosity > 3) 
          {
            printf("sample position: %g %g %g\n", pos.x, pos.y, pos.z); 
            printf("delta phi: %g\n",pp.source_phi - sum->peak[pol][peak].phi); 
            printf("delta el: %g\n",pp.source_theta - sum->peak[pol][peak].theta); 
            printf("payload el: %g %g\n", pp.payload_el, p.backwards_params.el_cutoff); 
          }
          //payload is below horizon.  or collides
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
        }
      

        /* compute the probabilities for each set of angles */ 
        pr.computeProbabilityDensity(nsamples, &phis[0], &thetas[0], &dens[0]); 

        /** Set occluded to 0
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
  //          printf("OCCLUDED\n"); 
            noccluded++; 
          }
          else if (dens[i] > max_dens)
          {
            max_dens = dens[i]; 
          }
        }

    //    printf("%d/%d samples occluded in segment %d\n", noccluded, nsamples, seg); 


        /* Now we want to compute the integral  */

        /** write out triangles for debugging if wanted*/ 
        if (debugfile)
        {
          if (!debugfile->Get("triangles")) debugfile->mkdir("triangles"); 
          debugfile->cd("triangles"); 
          TGraph2D g2d (nsamples,&phis[0], &thetas[0], &dens[0]); 
          g2d.Write(TString::Format("g%d",seg)); 
        }

  #ifdef HAVE_DELAUNAY
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

        double scale=  (max_phi-min_phi) * (max_theta-min_theta)  / ((del.XMax()-del.XMin())*(del.YMax()-del.YMin())); 

        double sum = 0; 
        bool enhance_flag = false;
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

        seg_p = sum; 
  //      printf("%d %g %d\n",seg, seg_p, noccluded); 
 

        if (p.verbosity > 2) 
        {

          printf("p(seg) %d: = %g\n",  seg, seg_p); 

        }

        if (seg_p > 1) 
        {
          if ( enhance_factor < p.backwards_params.max_enhance) 
          {
            enhance_flag = true; 
          }
          else
          {
            printf("OOPS %d %g %g\n", seg, scale, seg_p); 
            seg_p = 1; 
          }
        }

        done_with_this_segment = !enhance_flag; 
        enhance_factor++; //TODO, we can enhance more than once depending on how much above threshold we are 
      }

      if (occlusion) 
      {
          occlusion->push_back(std::pair<int,double> ( seg, double(noccluded)/nsamples)); 
      }



      contribution.push_back(std::pair<int,double>(seg,seg_p)); 
      if(max_densities) max_densities->push_back(std::pair<int,double>(seg, max_dens)); 
      /* if max density is above min_p, add the neighbors of this segment */ 
      if (max_dens >= min_p)
      {
        std::vector<int> new_neighbors;
        segmentationScheme()->getNeighbors(seg, &new_neighbors); 
        for (size_t j = 0; j < new_neighbors.size(); j++)
        {
          int new_seg = new_neighbors[j];
          if (!used[new_seg])
          {
            used[new_seg] = segs_to_check.size(); 
            segs_to_check.push_back(new_seg); 
          }
        }
      }
    }

  }

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
//         AntarcticCoord c(AntarcticCoord::WGS84,lat,lon,alt); 
         int seg = p.seg->getSegmentIndex(pp.source); 
         counts[seg]++; 
//          printf("%d %lld\n", seg, counts[seg]); 
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


  //finally, loop over any bases (if we were asked to do so) 
  if (base_contribution) 
  {

    //now loop over all the bases 
    int nbases = BaseList::getNumAbstractBases(); 
    for (int ibase = 0; ibase < nbases; ibase++)
    {

      const BaseList::abstract_base & base = BaseList::getAbstractBase(ibase); 
      AntarcticCoord base_pos= base.getPosition(gps->realTime);

      //TODO: can veto most bases early probalby 

      PayloadParameters pp (gps,base_pos, p.refract); 
      if (pp.payload_el < p.backwards_params.el_cutoff || ( p.collision_detection && pp.checkForCollision(p.collision_params.dx,0,0, p.dataset, p.collision_params.grace))) continue ; 

      double base_phi = pp.source_phi; 
      if (base_phi - phi0 > 180) base_phi-=360; 
      if (base_phi - phi0 < -180) base_phi+=360; 

      double dens = pr.computeProbabilityDensity( base_phi, pp.source_theta); 
 
      base_contribution->push_back(std::pair<int,double> (ibase, dens)); 
    }

  }
  
  return inv_two_pi_sqrt_det; 
}

int UCorrelator::ProbabilityMap::groupAdjacent(const double * vals_to_group, std::vector<std::vector<int> > * groups, double * counts, std::vector<double>  * dist) const
{
  std::vector<bool> consumed (p.seg->NSegments()); 
  int ngroups = 0; 


  //clean up input 
  if (groups) groups->clear(); 
  if (dist) dist->clear(); 
  if (counts) memset(counts, 0, sizeof(double) * p.seg->NSegments()); 

  for (int i = 0; i < p.seg->NSegments(); i++) 
  {
     if (! vals_to_group[i]) continue; 
     if (consumed[i]) continue; 

     ngroups++; 
     std::vector<int> group;
     group.push_back(i); 
     consumed[i] = true; 
     int n_in_group = 0; 
     double group_sum = vals_to_group[i]; 
     while (n_in_group < (int)  group.size()) 
     {
       int seg = group[n_in_group++]; 

       //see if we need to add any neighbors 

       std::vector<int> neighbors; 
       segmentationScheme()->getNeighbors(seg, &neighbors); 
       for (size_t i = 0; i < neighbors.size(); i++) 
       {
         if (vals_to_group[neighbors[i]] && !consumed[neighbors[i]])
         {
           consumed[neighbors[i]] = true;
           group.push_back(neighbors[i]); 
           group_sum += vals_to_group[neighbors[i]]; 
         }
       }
     }

     if (dist) dist->push_back(group_sum); 

     if (counts) 
     {
       for (size_t i = 0; i < group.size(); i++) 
       {
         counts[group[i]] = group_sum; 
       }
     }

     if (groups) groups->push_back(group); 
  }
  
  return ngroups; 

}


int UCorrelator::ProbabilityMap::makeMultiplicityTable(int level, bool blind, bool draw) const
{


  //non base table 
  std::vector<double> non_base; 

  std::vector<double> non_base_counts(draw ? segmentationScheme()->NSegments(): 0); 
  std::vector<double> wgt_without_base(getWgtAboveLevelWithoutBases(level), getWgtAboveLevelWithoutBases(level) + segmentationScheme()->NSegments()); 
  std::vector<double> wgt_with_base(getWgtAboveLevel(level), getWgtAboveLevel(level) + segmentationScheme()->NSegments()); 
  for (int i = 0; i < segmentationScheme()->NSegments(); i++) 
  {
    if (wgt_with_base[i] && wgt_without_base[i])
    {

      if (wgt_with_base[i] == wgt_without_base[i]) 
      {
        wgt_with_base[i] = 0; 
      }
      else
      {
        wgt_without_base[i] = 0; 
      }
    }
  }

  groupAdjacent(&wgt_without_base[0], 0,draw ? &non_base_counts[0] : 0, &non_base); 
  std::vector<double> base; 
  std::vector<double> base_counts(draw ? segmentationScheme()->NSegments(): 0); 
  groupAdjacent(&wgt_with_base[0], 0,draw ? &base_counts[0] : 0, &base); 

  const int maxes[] = {1,5,10,20,50,100,(int) 100e6}; 
  const int mins[] =  {1,2,6,11,21,51,101}; 
  int n_rows = sizeof(maxes) / sizeof(*maxes); 

  int n_base[n_rows] = { 0 }; 
  int n_non_base[n_rows] = { 0 } ; 


  for (size_t i = 0; i < non_base.size(); i++) 
  {
    for (int row = 0; row < n_rows; row++)
    {
       int n = round(non_base[i]); 
       if (n <= maxes[row] && n >= mins[row])
       {
         if (!(blind && row == 0))
         {
           n_non_base[row]++; 
         }
         break; 
       }
    }
  }

  for (size_t i = 0; i < base.size(); i++) 
  {
    for (int row = 0; row < n_rows; row++)
    {
       int n = round(base[i]); 
       if (n <= maxes[row] && n >= mins[row])
       {
         n_base[row]++; 
         break; 
       }
    }
  }




  printf("============================================================================\n"); 
  printf(" Clustering / base association based on Mahalanobis distance of %g\n", getLevel(level)); 
  if (blind) printf(" BLINDED TO NON-BASE SINGLES\n"); 

  printf("  Non-Base   |   Base  |   Size \n");
  printf("-------------------------------------\n");
  for (int row = 0; row < n_rows; row++) 
  {
    printf("   %04d    | %04d   |  ", n_non_base[row], n_base[row]);

    if (row == 0) 
      printf ("1\n"); 
    else if (row < n_rows-1)
      printf ("%d - %d\n", mins[row], maxes[row]); 
    else
      printf("%d - \n",mins[row]); 
  }
  printf("-------------------------------------\n"); 




  if (draw) 
  {
    for (int i = 0; i  < segmentationScheme()->NSegments(); i++) 
    {
      if (non_base_counts[i] && base_counts[i]) printf("OOPS!!!: %g %g\n", non_base_counts[i], base_counts[i]); 
      if (blind && non_base_counts[i] == 1)  non_base_counts[i] = 0; 
      non_base_counts[i] =  non_base_counts[i] - base_counts[i]; 
    }

    segmentationScheme()->Draw("colz",&non_base_counts[0]); 
  }

  return 0; 

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



