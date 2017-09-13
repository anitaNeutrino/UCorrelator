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



/** Convert a cumulative p value into instantaneous */ 
double UCorrelator::ProbabilityMap::cdf2density(double p) 
{
  return sqrt(-2 * log1p(-p)); // 
}

/** Convert a cumulative p value into instantaneous */ 
double UCorrelator::ProbabilityMap::density2cdf(double p) 
{
  return -expm1(-p*p/2.); 
}





static UCorrelator::ProbabilityMap::Params default_params; 

UCorrelator::ProbabilityMap::ProbabilityMap(const Params * par) 
  :  
    p( par ? *par : default_params), 
    ps(p.seg->NSegments(),0), 
    fraction_occluded(p.seg->NSegments(), 0), 
    NAboveLevel(p.n_level_thresholds, std::vector<int>(p.seg->NSegments(),0)),
    levels_p(p.n_level_thresholds),
    basesInSegment(p.seg->NSegments()), 
    baseNAboveLevel(p.n_level_thresholds, std::vector<int>(BaseList::getNumBases() + BaseList::getNumPaths())), 
    base_ps(BaseList::getNumBases() + BaseList::getNumPaths()) 
{

  for (size_t i = 0; i <NLevels(); i++)
  {
    levels_p[i] = cdf2density(p.level_cdf_thresholds[i]); 
  }

  for (size_t i = 0; i < BaseList::getNumBases(); i++)
  {
    const BaseList::base & b = BaseList::getBase(i); 
    int segment = p.seg->getSegmentIndex(b.getPosition(0));
    basesInSegment[segment].push_back(i); 
  }
}


int UCorrelator::ProbabilityMap::add(const AnitaEventSummary * sum, const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak, double weight, TFile * debugfile) 
{

  std::vector<std::pair<int,double> > segments_to_fill; 
  std::vector<std::pair<int,double> > base_ps_to_fill; 
  std::vector<std::pair<int,double> > occluded_to_fill; 

  computeContributions(sum,pat,pol,peak,segments_to_fill, &base_ps_to_fill, &occluded_to_fill, debugfile); 

  int NFilled = segments_to_fill.size(); 

  for (int i = 0; i < NFilled; i++)
  {
    int iseg = segments_to_fill[i].first; 
    double pseg = segments_to_fill[i].second; 
    ps[iseg] += pseg*weight; 
    fraction_occluded[iseg] += occluded_to_fill[i].second; 
    for (int j = 0; j < (int) NLevels(); j++)
    {
      if (pseg > levels_p[j]) 
      {
        NAboveLevel[j][iseg]++; 
      }
    }
  }

  int Nbases = base_ps_to_fill.size(); 
  for (int i = 0; i < Nbases; i++)
  {
    int ibase = base_ps_to_fill[i].first; 
    double pbase = base_ps_to_fill[i].second; 
    base_ps[ibase] += pbase; 

    for (int j = 0; j < (int) NLevels(); j++)
    {
      if (pbase > levels_p[j])
      {
        baseNAboveLevel[j][ibase]++; 
      }
    }
  }

  return NFilled; 
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


  if (minDensity() != other.minDensity())
  {
    fprintf(stderr,"Combining ProbabilityMap's with different min densities (ours=%g,theirs=%g). Continuing, but may not be what you want!\n", minDensity(), other.minDensity());  
  }
  

  //TODO we probably want to compare the pointing resolution models as well 

  //if we made it this far, all is good (except for what hasn't been implemented yet) !

  for (int i = 0; i < segmentationScheme()->NSegments(); i++) 
  {
    ps[i] += other.getDensitySums()[i]; 

    for (int j = 0; j < (int) NLevels(); j++)
    {
      NAboveLevel[j][i] += other.getNAboveLevel(j)[i]; 
    }
  }

  if (combine_bases) 
  {
    for (size_t i = 0; i < getNBases(); i++) 
    {

      base_ps[i] += other.getBaseDensitySums()[i]; 
      for (int j = 0; j < (int) NLevels(); j++)
      {
        baseNAboveLevel[j][i] += other.getBaseNAboveLevel(j)[i]; 
      }
    }
  }

  return 0; 
}


void UCorrelator::ProbabilityMap::computeContributions(const AnitaEventSummary * sum, const Adu5Pat * gps, 
                                                       AnitaPol::AnitaPol_t pol, int peak,
                                                      std::vector<std::pair<int,double> > & contribution, 
                                                      std::vector<std::pair<int,double> > * base_contribution, 
                                                      std::vector<std::pair<int,double> > * occlusion, 
                                                      TFile * debugfile) const 
{

  /* compute pointing resolution */ 
  PointingResolution pr; 
  p.point->computePointingResolution(sum,pol, peak, &pr);  
  contribution.clear(); 

  std::vector<int> used ( segmentationScheme()->NSegments()); 
  UsefulAdu5Pat pat(gps); 
  if (p.projection == Params::BACKWARD) 
  {
    //start with guess
    const AnitaEventSummary::PointingHypothesis *pk = &sum->peak[pol][peak];  
    AntarcticCoord guess(AntarcticCoord::WGS84, pk->latitude, pk->longitude, pk->altitude); 

    PayloadParameters check(gps,guess); 

    printf("payload coords: %g %g\n", sum->peak[pol][peak].phi, sum->peak[pol][peak].theta); 
    printf("guess coords: %g %g %g %g %g\n", check.source_phi, check.source_theta, check.payload_az, check.payload_el, check.distance); 
    printf("guess position: %g %g %g\n", guess.x, guess.y, guess.z); 
   
    /* set up vector of segments to check */
    size_t nchecked = 0; 
    std::vector<int> segs_to_check; 
    segs_to_check.reserve(100);  // a plausible number
    int guess_seg = segmentationScheme()->getSegmentIndex(guess); 
    if (guess_seg < 0) return; 
    used[guess_seg] = nchecked; 
    segs_to_check.push_back(guess_seg); 

    AntarcticCoord pos = segmentationScheme()->getSegmentCenter(segs_to_check[0]); 
    pos.to(AntarcticCoord::WGS84); 
    printf("center position: %g %g %g\n", pos.x, pos.y, pos.z); 

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


      //check for any stationary bases here 
      int nbases = basesInSegment[seg].size(); 
      if (base_contribution && nbases > 0)
      {

        //first, handle stationary bases 

        std::vector<double> base_phis(nbases); 
        std::vector<double> base_thetas(nbases); 
        std::vector<double> base_dens(nbases); 
        std::vector<int> base_list(nbases); 
        //loop over the bases within this segment 
        int nbases_used = 0;
        for (int ibase = 0; ibase < nbases; ibase++)
        {
          int basenum = basesInSegment[seg][ibase]; 
          const BaseList::base & base =  BaseList::getBase(basenum); 
          printf("Considering base %d: %s\n",basenum, base.getName()); 

          PayloadParameters pp (gps,base.getPosition(gps->realTime)); 

          if (pp.payload_el < p.backwards_params.el_cutoff || (p.collision_detection && pp.checkForCollision(p.collision_params.dx,0, p.dataset, p.collision_params.grace)))
          {
            continue; 
          }

          base_phis[nbases_used] = pp.source_phi; 
          base_thetas[nbases_used] = pp.source_theta; 
          base_list[nbases_used] = basenum; 
          nbases_used++; 
        }

        pr.computeProbabilityDensity(nbases_used, &base_phis[0], &base_thetas[0], &base_dens[0]); 

        for (int i = 0; i < nbases_used; i++)
        {
          base_contribution->push_back(std::pair<int,double> (base_list[i], base_dens[i])); 
        }


      }
      


  //    printf("%g %g %g %g %g\n", check.source_phi, check.source_theta, check.payload_az, check.payload_el, check.distance); 
      
      // done with bases, now we sample our 
      // segment into a bunch of positions


      segmentationScheme()->sampleSegment(seg, nsamples, &samples[0], p.backwards_params.random_samples);  // do we want to randomize? i dunno. I'll decide later. 
      
      // loop over the samples, 
      //  we want to check:
      //     - can this sample probably see ANITA? 
      //     - what are the coordinates of this sample in ANITA's frame? 
      for (int i = 0; i < nsamples; i++)
      {
        //This computes the payload in source coords and vice versa
        PayloadParameters pp(gps, samples[i]); 

        pos = samples[i].as(AntarcticCoord::STEREOGRAPHIC); 
//        printf("sample position: %g %g %g\n", pos.x, pos.y, pos.z); 
//        printf("delta phi: %g\n",pp.source_phi - sum->peak[pol][peak].phi); 
//        printf("delta el: %g\n",pp.source_theta - sum->peak[pol][peak].theta); 
//        printf("payload el: %g\n", pp.payload_el); 
  //
        //payload is below horizon. 

        AntarcticCoord collid; 
        if (pp.payload_el < p.backwards_params.el_cutoff || ( p.collision_detection && pp.checkForCollision(p.collision_params.dx,&collid, p.dataset, p.collision_params.grace))) 
        {
          if (pp.payload_el >=p.backwards_params.el_cutoff) 
          {
 //           collid.to(AntarcticCoord::STEREOGRAPHIC); 
//            printf("Collision at (%g %g %g), surface at %g\n", collid.x,collid.y,collid.z, RampdemReader::SurfaceAboveGeoidEN(collid.x,collid.y, p.dataset)); 


             // we want to make sure that whatever segment is occluding is is considered, if it isn't already. So let's project to continent from payload and ensure we have that segment already 
             AntarcticCoord project_to(AntarcticCoord::WGS84,0,0,0); 
             int success = pat.getSourceLonAndLatAtAlt(pp.source_phi * TMath::DegToRad(), pp.source_theta * TMath::DegToRad(), project_to.x, project_to.y, project_to.z); 
             if (success == 1) 
             {
               int potential_seg = p.seg->getSegmentIndex(project_to); 
               if (!used[potential_seg]) 
               {
                 used[potential_seg] = segs_to_check.size(); 
                 segs_to_check.push_back(potential_seg); 
               }
             }
            
          }
          occluded[i] = true; 
        }

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
      for (int i =0; i < nsamples; i++)
      {
        if (occluded[i])
        {
          dens[i] = 0; 
          noccluded++; 
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
  //    printf("%d %g\n",seg, seg_p); 

      contribution.push_back(std::pair<int,double>(seg,seg_p)); 
      /* if p is above min_p, add the neighbors of this segment */ 
      if (seg_p >= minDensity())
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
      double lon,lat,alt; 
//        printf("%f %f\n",phi,theta); 
      int success = pat.getSourceLonAndLatAtAlt(phi * TMath::DegToRad(), theta * TMath::DegToRad(), lon,lat,alt); 

      if (success == 1) 
      {
         AntarcticCoord c(AntarcticCoord::WGS84,lat,lon,alt); 
         int seg = p.seg->getSegmentIndex(c); 
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

    //now loop over all the bases 
    int nbases = BaseList::getNumBases(); 
    for (int ibase = 0; ibase < nbases; ibase++)
    {

      const BaseList::base & base = BaseList::getBase(ibase); 
      PayloadParameters pp (gps,base.getPosition(gps->realTime)); 
      if (p.collision_detection && pp.checkForCollision(p.collision_params.dx,0, p.dataset, p.collision_params.grace))  continue; 
      double dens = pr.computeProbabilityDensity( pp.source_phi, pp.source_theta); 
      base_contribution->push_back(std::pair<int,double> (ibase, dens)); 
    }
  }


  //finally, loop over any paths that are active 

  if (base_contribution) 
  {

    int npaths = BaseList::getNumPaths(); 
    int path_offset = BaseList::getNumBases(); 
    for (int ipath = 0; ipath < npaths; ipath++)
    {
      const BaseList::path & path = BaseList::getPath(ipath); 
      if (!path.isValid(gps->realTime)) continue; //skip invalid paths at time

      PayloadParameters pp (gps,path.getPosition(gps->realTime)); 

      //because the path may be above the ground, we always have to do a collision check. 
      if (p.collision_detection && pp.checkForCollision(p.collision_params.dx,0, p.dataset, p.collision_params.grace))  continue; 
      double dens = pr.computeProbabilityDensity( pp.source_phi, pp.source_theta); 
      base_contribution->push_back(std::pair<int,double> (path_offset + ipath, dens)); 
    }
  }
  
}



