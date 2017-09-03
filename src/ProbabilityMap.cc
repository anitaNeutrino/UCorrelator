#include "ProbabilityMap.h" 
#include "PointingResolutionModel.h"
#include "Adu5Pat.h" 
#include "BaseList.h" 
 

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
  return sqrt(-2 * log(1-p)); 
}

/** Convert a cumulative p value into instantaneous */ 
double UCorrelator::ProbabilityMap::density2cdf(double p) 
{
  return 1 - exp(-p*p/2.); 
}





static StereographicGrid defaultSeg; 
static UCorrelator::ConstantPointingResolutionModel defaultPoint; 

UCorrelator::ProbabilityMap::ProbabilityMap(const AntarcticSegmentationScheme * seg, 
                     const PointingResolutionModel * pointing_model, 
                     int NLevelThresholds, const double * level_thresholds, 
                     double cutoff, int num_samples_per_bin) 
  :  g(seg ? seg : &defaultSeg), 
    prm(pointing_model ? pointing_model : &defaultPoint),
    ps(g->NSegments(),0), 
    NAboveLevel(NLevelThresholds, std::vector<int>(g->NSegments(),0)),
    levels(NLevelThresholds), levels_p(NLevelThresholds), 
    min_p(cdf2density(cutoff)),cutoff(cutoff), 
    basesInSegment(g->NSegments()), 
    baseNAboveLevel(NLevelThresholds, std::vector<int>(BaseList::getNumBases() + BaseList::getNumPaths())), 
    base_ps(BaseList::getNumBases() + BaseList::getNumPaths()), 
    nsamples(num_samples_per_bin)
    
{
  memcpy(&levels[0], level_thresholds, NLevelThresholds*sizeof(double)); 

  for (size_t i = 0; i < levels.size(); i++)
  {
    levels_p[i] = cdf2density(levels[i]); 
  }

  for (size_t i = 0; i < BaseList::getNumBases(); i++)
  {
    const BaseList::base & b = BaseList::getBase(i); 
    int segment = g->getSegmentIndex(b.getPosition(0));
    basesInSegment[segment].push_back(i); 
  }
}



int UCorrelator::ProbabilityMap::add(const AnitaEventSummary * sum, const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak, double weight, TFile * debugfile) 
{

  std::vector<std::pair<int,double> > segments_to_fill; 
  std::vector<std::pair<int,double> > base_ps_to_fill; 

  computeContributions(sum,pat,pol,peak,segments_to_fill, &base_ps_to_fill, debugfile); 

  int NFilled = segments_to_fill.size(); 
  for (int i = 0; i < NFilled; i++)
  {
    int iseg = segments_to_fill[i].first; 
    double pseg = segments_to_fill[i].second; 
    ps[iseg] += pseg*weight; 
    for (int j = 0; j < (int) levels.size(); j++)
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

    for (int j = 0; j < (int) levels.size(); j++)
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

  g->asString(&our_scheme); 
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


  if (getCutoff() != other.getCutoff())
  {
    fprintf(stderr,"Combining ProbabilityMap's with different cutoffs (ours=%g,theirs=%g). Continuing, but may not be what you want!\n", getCutoff(), other.getCutoff());  
  }
  

  //TODO we probably want to compare the pointing resolution models as well 

  //if we made it this far, all is good (except for what hasn't been implemented yet) !

  for (int i = 0; i < g->NSegments(); i++) 
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
                                                      TFile * debugfile) const 
{

  //start with guess
  const AnitaEventSummary::PointingHypothesis *pk = &sum->peak[pol][peak];  
  AntarcticCoord guess(AntarcticCoord::WGS84, pk->latitude, pk->longitude, pk->altitude); 

  PayloadParameters check(gps,guess); 
//  printf("payload coords: %g %g\n", sum->peak[pol][peak].phi, sum->peak[pol][peak].theta); 
//  printf("guess coords: %g %g %g %g %g\n", check.source_phi, check.source_theta, check.payload_az, check.payload_el, check.distance); 
//  printf("guess position: %g %g %g\n", guess.x, guess.y, guess.z); 
 
  /* set up lvector of segments to check */
  size_t nchecked = 0; 
  std::vector<int> segs_to_check; 
  std::vector<bool> used(g->NSegments()); 
  segs_to_check.reserve(100);  // a plausible number
  int guess_seg = g->getSegmentIndex(guess); 
  if (guess_seg < 0) return; 
  segs_to_check.push_back(guess_seg); 
  used[guess_seg] = true; 

  AntarcticCoord pos = g->getSegmentCenter(segs_to_check[0]); 
  pos.to(AntarcticCoord::WGS84); 
//  printf("center position: %g %g %g\n", pos.x, pos.y, pos.z); 
   


  /* compute pointing resolution */ 
  PointingResolution p; 
  prm->computePointingResolution(sum,pol, peak, &p);  

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
      std::vector<int> base_list; 
      //loop over the bases within this segment 
      int nbases_used = 0;
      for (int ibase = 0; ibase < nbases; ibase++)
      {
        int basenum = basesInSegment[seg][ibase]; 
        const BaseList::base & base =  BaseList::getBase(basenum); 

        //TODO for transient bases, check if it's active at this point, once the base has a way of letting us know
        //
        PayloadParameters pp (gps,base.getPosition(gps->realTime)); 

        if (pp.payload_el < 0 || pp.checkForCollision())
        {
          continue; 
        }

        base_phis[nbases_used] = pp.source_phi; 
        base_thetas[nbases_used] = pp.source_theta; 
        base_list[nbases_used] = basenum; 
        nbases_used++; 
      }

      p.computeProbabilityDensity(nbases_used, &base_phis[0], &base_thetas[0], &base_dens[0]); 

      for (int i = 0; i < nbases_used; i++)
      {
        base_contribution->push_back(std::pair<int,double> (base_list[i], base_dens[i])); 
      }


    }
    


//    printf("%g %g %g %g %g\n", check.source_phi, check.source_theta, check.payload_az, check.payload_el, check.distance); 
    
    // done with bases, now we sample our 
    // segment into a bunch of positions

    g->sampleSegment(seg, nsamples, &samples[0],false);  // do we want to randomize? i dunno. I'll decide later. 
    
    // loop over the samples, 
    //  we want to check:
    //     - can this sample probably see ANITA? 
    //     - what are the coordinates of this sample in ANITA's frame? 
    for (int i = 0; i < nsamples; i++)
    {
      //This computes the payload in source coords and vice versa
      PayloadParameters pp(gps, samples[i]); 

 //     pos = samples[i].as(AntarcticCoord::WGS84); 
//      printf("sample position: %g %g %g\n", pos.x, pos.y, pos.z); 
//      printf("%g\n",pp.source_phi - sum->peak[pol][peak].phi); 
//
      //payload is below horizon. I should actually compare to local gradient instead of 0... 
      AntarcticCoord coll;
      if (pp.payload_el < 0 || pp.checkForCollision(100,&coll)) 
      {
        TString str;
        coll.as(AntarcticCoord::STEREOGRAPHIC).asString(&str); 
        printf("Collision at %s\n", str.Data()); 
        occluded[i] = true; 
      }


      /* add to phis and thetas for this segment */
      phis[i]=pp.source_phi; 
      thetas[i]=pp.source_theta; 
    }
    

    /* compute the probabilities for each set of angles */ 
    p.computeProbabilityDensity(nsamples, &phis[0], &thetas[0], &dens[0]); 

    /** Set occluded to 0
     *
     * Why not just get rid of them entirely you ask? I think the zero's might be important
     * for calculating the integral properly. But maybe they're not. I'll look into it. 
     *
     * */ 
    for (int i =0; i < nsamples; i++)
    {
      if (occluded[i]) dens[i] = 0; 
    }

    /* Now we want to integral  */


    /** write out triangles for debugging if wanted*/ 
    if (debugfile)
    {
      debugfile->cd(); 
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

    /* if p is above min_p, include this contribution and add the neighbors of this segment */ 
    if (seg_p > min_p)
    {
      contribution.push_back(std::pair<int,double>(seg,seg_p)); 
      std::vector<int> new_neighbors;
      g->getNeighbors(seg, &new_neighbors); 
      for (size_t j = 0; j < new_neighbors.size(); j++)
      {
        int new_seg = new_neighbors[j];
        if (!used[new_seg])
        {
          segs_to_check.push_back(new_seg); 
          used[new_seg] = true; 
        }
      }
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
      if (pp.checkForCollision())  continue; 

      double dens = p.computeProbabilityDensity( pp.source_phi, pp.source_theta); 
      base_contribution->push_back(std::pair<int,double> (path_offset + ipath, dens)); 
    }
  }
  
}



