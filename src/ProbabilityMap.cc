#include "ProbabilityMap.h" 
#include "PointingResolutionModel.h"
#include "Adu5Pat.h" 
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
    ps_with_base(p.seg->NSegments(),0), 
    ps_without_base(p.seg->NSegments(),0), 
    ps_norm(p.seg->NSegments(),0), 
    ps_norm_with_base(p.seg->NSegments(),0), 
    ps_norm_without_base(p.seg->NSegments(),0), 
    fraction_occluded(p.seg->NSegments(), 0), 
    n_above_level(p.n_level_thresholds, std::vector<int>(p.seg->NSegments(),0)),
    wgt_above_level(p.n_level_thresholds, std::vector<double>(p.seg->NSegments(),0)),
    bases_in_segment(p.seg->NSegments()), 
    n_above_level_with_base(p.n_level_thresholds, std::vector<int>(p.seg->NSegments(),0)),
    wgt_above_level_with_base(p.n_level_thresholds, std::vector<double>(p.seg->NSegments(),0)),
    n_above_level_without_base(p.n_level_thresholds, std::vector<int>(p.seg->NSegments(),0)),
    wgt_above_level_without_base(p.n_level_thresholds, std::vector<double>(p.seg->NSegments(),0)),
    base_n_above_level(p.n_level_thresholds, std::vector<int>(BaseList::getNumBases() + BaseList::getNumPaths(), 0)), 
    base_ps(BaseList::getNumBases() + BaseList::getNumPaths(), 0.) 
{


  for (size_t i = 0; i < BaseList::getNumBases(); i++)
  {
    const BaseList::base & b = BaseList::getBase(i); 
    int segment = p.seg->getSegmentIndex(b.getPosition(0));
    if (segment !=-1) 
      bases_in_segment[segment].push_back(i); 
  }
}


int UCorrelator::ProbabilityMap::add(const AnitaEventSummary * sum, const Adu5Pat * pat, AnitaPol::AnitaPol_t pol,
                                     int peak, double weight, TFile * debugfile) 
{

  std::vector<std::pair<int,double> > segments_to_fill; 
  std::vector<std::pair<int,double> > base_ps_to_fill; 
  std::vector<std::pair<int,double> > occluded_to_fill; 
  std::vector<std::pair<int,double> > max_densities; ; 

  double inv_two_pi_sqrt_det = computeContributions(sum,pat,pol,peak,segments_to_fill, &base_ps_to_fill, &occluded_to_fill, &max_densities, debugfile); 

  std::vector<double> levels_p(p.n_level_thresholds); 

  for (int i = 0; i < (p.n_level_thresholds); i++)
  {
    levels_p[i] = dist2dens(p.level_thresholds[i], inv_two_pi_sqrt_det); 
//    printf("mahalanobis distance: %g, probability density: %g\n", p.level_thresholds[i], levels_p[i]); 
  }

  int incr = weight > 0 ? 1 : weight < 0 ? -1 : 0; 
  int Nbases = base_ps_to_fill.size(); 
  double total_base_level = 0; 
  for (int i = 0; i < Nbases; i++)
  {
    int ibase = base_ps_to_fill[i].first; 
    double pbase = base_ps_to_fill[i].second; 
    base_ps[ibase] += pbase * weight; 
    total_base_level += pbase * weight; 
    printf("Base contribution from %s (%d) : %g, total is %g:\n", BaseList::getAbstractBase(ibase).getName(), ibase, pbase,  base_ps[ibase]); 

    for (int j = 0; j < (int) NLevels(); j++)
    {
      if (pbase > levels_p[j])
      {
        base_n_above_level[j][ibase]+=incr; 
      }
    }

//    printf("Base contribution from %s (%d) : %g, total is %g:\n", BaseList::getAbstractBase(ibase).getName(), ibase, pbase,  base_ps[ibase]); 
  }



  int NFilled = segments_to_fill.size(); 


  double total_ps = 0; 
  std::vector<int> this_NAboveLevel(NLevels()); 

  for (int i = 0; i < NFilled; i++)
  {
    double pseg = segments_to_fill[i].second; 
    total_ps += pseg; 
    double max_dens_seg = max_densities[i].second; 
    for (int j = 0; j < (int) NLevels(); j++) 
    {
      if (max_dens_seg > levels_p[j])
      {
        this_NAboveLevel[j]+=incr; 
      }
    }
  }


  for (int i = 0; i < NFilled; i++)
  {
    int iseg = segments_to_fill[i].first; 
    double pseg = segments_to_fill[i].second; 
    double max_dens_seg = max_densities[i].second; 
    ps[iseg] += pseg*weight; 
    ps_norm[iseg] += pseg * weight / total_ps; 
    if (total_base_level > 0) 
    {
      ps_with_base[iseg] += pseg*weight; 
      ps_norm_with_base[iseg] += pseg * weight / total_ps; 
    }
    else
    {
      ps_without_base[iseg] += pseg*weight; 
      ps_norm_without_base[iseg] += pseg * weight / total_ps; 
    }

    fraction_occluded[iseg] += occluded_to_fill[i].second; 
    for (int j = 0; j < (int) NLevels(); j++)
    {
      if (max_dens_seg > levels_p[j]) 
      {
//        printf("Segment %d above level %d (%g) \n", iseg, j, pseg); 
        n_above_level[j][iseg]+=incr; 
        wgt_above_level[j][iseg]+= incr/ double(this_NAboveLevel[j]); 
        if (total_base_level > levels_p[j])
        {
          n_above_level_with_base[j][iseg]+=incr; 
          wgt_above_level_with_base[j][iseg]+= incr/ double(this_NAboveLevel[j]); 
        }
        else
        {
          n_above_level_without_base[j][iseg]+=incr; 
          wgt_above_level_without_base[j][iseg]+= incr/ double(this_NAboveLevel[j]); 
        }
      }
    }
  }

//  dumpNonZeroBases(); 



  return NFilled; 
}

int UCorrelator::ProbabilityMap::dumpNonZeroBases()  const 
{
  int n = 0;
  for (int i = 0; i < (int) getNBases(); i++)
  {
    if(base_ps[i] > 0) 
    {
      printf("%s (%d): %g\n", BaseList::getAbstractBase(i).getName(), i, base_ps[i]); 
      n++; 
    }
  }
 
  return n; 
}

double UCorrelator::ProbabilityMap::overlap(const AnitaEventSummary * sum , const Adu5Pat * pat,  AnitaPol::AnitaPol_t pol, int peak, 
                                            std::vector<std::pair<int,double> > * bases, bool remove_self) const
{
  std::vector<std::pair<int,double> > segs; 

  computeContributions(sum,pat,pol,peak,segs,bases); 

  int N = segs.size(); 
  double answer = 0;

  for (int i =0; i< N; i++)
  {
    answer += ps[segs[i].first]; 
    if (remove_self) answer -= segs[i].second; 
  }


  
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
    ps[i] += other.getDensitySums()[i]; 
    ps_norm[i] += other.getDensitySumsNormalized()[i]; 
    ps_with_base[i] += other.getDensitySumsWithBases()[i]; 
    ps_norm_with_base[i] += other.getDensitySumsNormalizedWithBases()[i]; 
    ps_without_base[i] += other.getDensitySumsWithoutBases()[i]; 
    ps_norm_without_base[i] += other.getDensitySumsNormalizedWithoutBases()[i]; 

    for (int j = 0; j < (int) NLevels(); j++)
    {
      n_above_level[j][i] += other.getNAboveLevel(j)[i]; 
      wgt_above_level[j][i] += other.getWgtAboveLevel(j)[i]; 
      n_above_level_with_base[j][i] += other.getNAboveLevelWithBases(j)[i]; 
      wgt_above_level_with_base[j][i] += other.getWgtAboveLevelWithBases(j)[i]; 
      n_above_level_without_base[j][i] += other.getNAboveLevelWithoutBases(j)[i]; 
      wgt_above_level_without_base[j][i] += other.getWgtAboveLevelWithoutBases(j)[i]; 
    }
  }

  if (combine_bases) 
  {
    for (size_t i = 0; i < getNBases(); i++) 
    {

      base_ps[i] += other.getBaseDensitySums()[i]; 
      for (int j = 0; j < (int) NLevels(); j++)
      {
        base_n_above_level[j][i] += other.getBaseNAboveLevel(j)[i]; 
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

  double inv_two_pi_sqrt_det = get_inv_two_pi_sqrt_det(pr.getdPhi(), pr.getdTheta(), pr.getCorr()); 
  double min_p = dist2dens(maxDistance(), inv_two_pi_sqrt_det); 

  std::vector<int> used ( segmentationScheme()->NSegments()); 
  UsefulAdu5Pat pat(gps); 

  if (p.projection == Params::BACKWARD) 
  {
    //start with guess
    const AnitaEventSummary::PointingHypothesis *pk = &sum->peak[pol][peak];  
    PayloadParameters guess;  
    PayloadParameters::findSourceOnContinent(pk->theta,pk->phi,gps, &guess, p.refract, p.collision_detection ? p.collision_params.dx : 0); 
//    guess.source.to(AntarcticCoord::STEREOGRAPHIC); 
//    printf("%d %g %g %g\n", status, guess.source.x, guess.source.y, guess.source.z); 
   
    /* set up vector of segments to check */
    size_t nchecked = 0; 
    std::vector<int> segs_to_check; 
    segs_to_check.reserve(100);  // a plausible number
    int guess_seg = segmentationScheme()->getSegmentIndex(guess.source); 
    if (guess_seg < 0) return inv_two_pi_sqrt_det; 
    used[guess_seg] = nchecked; 
    segs_to_check.push_back(guess_seg); 

    AntarcticCoord pos = segmentationScheme()->getSegmentCenter(segs_to_check[0]); 
    pos.to(AntarcticCoord::WGS84); 
 //   printf("center position: %g %g %g\n", pos.x, pos.y, pos.z); 

    int nsamples = p.backwards_params.num_samples_per_bin;
    /* set up vectors for sample phis /thetas/ densities */ 
    std::vector<AntarcticCoord> samples(nsamples); 
    std::vector<double> dens(nsamples);; 
    std::vector<double> phis(nsamples);; 
    std::vector<double> thetas(nsamples);; 
    std::vector<bool> occluded(nsamples); 

    /** Loop over segments we need to check to see if p > cutoff */
    while (nchecked < segs_to_check.size())
    {
      int seg = segs_to_check[nchecked++]; 


  //    printf("%g %g %g %g %g\n", check.source_phi, check.source_theta, check.payload_az, check.payload_el, check.distance); 
      
      // segment into a bunch of positions
      segmentationScheme()->sampleSegment(seg, nsamples, &samples[0], p.backwards_params.random_samples);  // do we want to randomize? i dunno. I'll decide later. 
      
      // loop over the samples, 
      //  we want to check:
      //     - can this sample probably see ANITA? 
      //     - what are the coordinates of this sample in ANITA's frame? 
      for (int i = 0; i < nsamples; i++)
      {
        //This computes the payload in source coords and vice versa
        PayloadParameters pp(gps, samples[i], p.refract); 

        pos = samples[i].as(AntarcticCoord::STEREOGRAPHIC); 
//        printf("sample position: %g %g %g\n", pos.x, pos.y, pos.z); 
//        printf("delta phi: %g\n",pp.source_phi - sum->peak[pol][peak].phi); 
//        printf("delta el: %g\n",pp.source_theta - sum->peak[pol][peak].theta); 
//        printf("payload el: %g %g\n", pp.payload_el, p.backwards_params.el_cutoff); 
  //
        //payload is below horizon.  or collides
        AntarcticCoord collid, collid_exit; 
        if (pp.payload_el < p.backwards_params.el_cutoff || ( p.collision_detection && pp.checkForCollision(p.collision_params.dx,&collid,&collid_exit, p.dataset, p.collision_params.grace))) 
        {
          if (pp.payload_el >=p.backwards_params.el_cutoff)  //colides
          {
             collid.to(AntarcticCoord::STEREOGRAPHIC); 
//             printf("Collision at (%g %g %g), surface at %g\n", collid.x,collid.y,collid.z, RampdemReader::SurfaceAboveGeoidEN(collid.x,collid.y, p.dataset)); 


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
      int noccluded = 0; 
      double max_dens = 0;
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

      if (occlusion) 
      {
        occlusion->push_back(std::pair<int,double> ( seg, double(noccluded)/nsamples)); 
      }

  //    printf("%d/%d samples occluded in segment %d\n", noccluded, nsamples, seg); 


      /* Now we want to compute the integral  */

      /** write out triangles for debugging if wanted*/ 
      if (debugfile)
      {
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
      for (std::vector<ROOT::Math::Delaunay2D::Triangle>::const_iterator it = del.begin(); it!=del.end(); it++)
      {
        const ROOT::Math::Delaunay2D::Triangle & tri = *it; 

        double area = 0.5 *scale* fabs(( tri.x[0] - tri.x[2]) * (tri.y[1] - tri.y[0]) - (tri.x[0] - tri.x[1]) * (tri.y[2] - tri.y[0])); 
        for (int j = 0; j < 3; j++) 
        {
          sum += area /3. *  dens[tri.idx[j]];
        }
      }
#else
      double sum =-1; 
      fprintf(stderr,"ROOT 5 not currently supported in Probability Map due to lack of ROOT/Math/Delaunay2D.h. Will be fixed someday if necessary. \n"); 
#endif

      double seg_p = sum; 
//      printf("%d %g %d\n",seg, seg_p, noccluded); 

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
      PayloadParameters pp (gps,base.getPosition(gps->realTime), p.refract); 
      if (pp.payload_el < p.backwards_params.el_cutoff || ( p.collision_detection && pp.checkForCollision(p.collision_params.dx,0,0, p.dataset, p.collision_params.grace)))  continue; 
      double dens = pr.computeProbabilityDensity( pp.source_phi, pp.source_theta); 
      base_contribution->push_back(std::pair<int,double> (ibase, dens)); 
    }

  }
  
  return inv_two_pi_sqrt_det; 
}

int UCorrelator::ProbabilityMap::groupAdjacent(const double * vals_to_group, double * counts, std::vector<double>  * dist) const
{
  std::vector<bool> consumed (p.seg->NSegments()); 
  int ngroups = 0; 
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

     for (size_t i = 0; i < group.size(); i++) counts[group[i]] = group_sum; 
  }
  
  return ngroups; 

}




