#include "ProbabilityMap.h" 
#include "PointingResolutionModel.h"
#include "Adu5Pat.h" 

ClassImp(UCorrelator::ProbabilityMap); 

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
    levels(NLevelThresholds), min_p(cutoff), nsamples(num_samples_per_bin)
{
  memcpy(&levels[0], level_thresholds, NLevelThresholds*sizeof(double)); 
}





int UCorrelator::ProbabilityMap::add(const AnitaEventSummary * sum, const Adu5Pat * pat, AnitaPol::AnitaPol_t pol, int peak, double weight) 
{

  std::vector<std::pair<int,double> > segments_to_fill; 
  computeContributions(sum,pat,pol,peak,segments_to_fill); 

  int NFilled = segments_to_fill.size(); 
  for (int i = 0; i < NFilled; i++)
  {
    int iseg = segments_to_fill[i].first; 
    double pseg = segments_to_fill[i].second; 
    ps[iseg] += pseg*weight; 
    for (int j = 0; j < (int) levels.size(); j++)
    {
      if (pseg > levels[j]) 
      {
        NAboveLevel[j][i]++; 
      }
    }
  }

  return NFilled; 
}


double UCorrelator::ProbabilityMap::overlap(const AnitaEventSummary * sum , const Adu5Pat * pat,  AnitaPol::AnitaPol_t pol, int peak, bool remove_self) const
{
  std::vector<std::pair<int,double> > segs; 
  computeContributions(sum,pat,pol,peak,segs); 

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


  if (probabilityCutoff() != other.probabilityCutoff())
  {
    fprintf(stderr,"Combining ProbabilityMap's with different cutoffs (ours=%g,theirs=%g). Continuing, but may not be what you want!\n", probabilityCutoff(), other.probabilityCutoff());  
  }
  
  //TODO we probably want to compare the pointing resolution models as well 

  //if we made it this far, all is good (except for what hasn't been implemented yet) !

  for (int i = 0; i < g->NSegments(); i++) 
  {
    ps[i] += other.getProbabilities()[i]; 

    for (int j = 0; j < (int) NLevels(); j++)
    {
      NAboveLevel[j][i] += other.getNAboveLevel(j)[i]; 
    }
  }


  return 0; 
}


void UCorrelator::ProbabilityMap::computeContributions(const AnitaEventSummary * sum, const Adu5Pat * gps, AnitaPol::AnitaPol_t pol, int peak,  std::vector<std::pair<int,double> > & contribution) const 
{

  //start with guess
  const AnitaEventSummary::PointingHypothesis *pk = &sum->peak[pol][peak];  
  AntarcticCoord guess(AntarcticCoord::WGS84, pk->latitude, pk->longitude, pk->altitude); 
 
  /* set up lvector of segments to check */
  size_t nchecked = 0; 
  std::vector<int> segs_to_check; 
  segs_to_check.reserve(100);  // a plausible number
  segs_to_check.push_back(g->getSegmentIndex(guess)); 

  /* compute pointing resolution */ 
  PointingResolution p; 
  prm->computePointingResolution(sum,pol, peak, &p);  

  /* set up vectors for sample phis /thetas/ probs */ 
  std::vector<AntarcticCoord> samples(nsamples); 
  std::vector<double> probs; 
  std::vector<double> phis; 
  std::vector<double> thetas; 
  probs.reserve(nsamples); 
  phis.reserve(nsamples); 
  thetas.reserve(nsamples); 

  /** Loop over segments we need to check to see if p > cutoff */
  while (nchecked < segs_to_check.size())
  {
    int seg = segs_to_check[nchecked++]; 

    probs.clear(); 
    phis.clear(); 
    thetas.clear(); 

    // sample our segment into a bunch of positions
    g->sampleSegment(seg, nsamples, &samples[0],false);  // do we want to randomize? i dunno. I'll decide later. 
    

    // loop over the samples, 
    //  we want to check:
    //     - can this sample probably see ANITA? 
    //     - what are the coordinates of this sample in ANITA's frame? 
    for (int i = 0; i < nsamples; i++)
    {
      //This computes the payload in source coords and vice versa
      PayloadParameters pp(gps, samples[i]); 

      //payload is below horizon. I should actually compare to local gradient, probably, or do some ray tracing here
      if (pp.payload_el < 0)
         continue; 


      /* add to phis and thetas for this segment */
      phis.push_back(pp.source_phi); 
      thetas.push_back(pp.source_theta); 
    }
    
    /* compute mean of probabilities for this segment. Is this correct? I'm not sure */ 
    probs.resize(phis.size()); 
    p.computeProbability(phis.size(), &phis[0], &thetas[0], &probs[0]); 
    double seg_p = TMath::Mean(probs.size(), &probs[0]); 

    /* if p is above cutoff, include this contribution and add the neighbors of this segment */ 
    if (seg_p > min_p) 
    {
      contribution.push_back(std::pair<int,double>(seg,seg_p)); 
      g->getNeighbors(seg, &segs_to_check); 
    }
  }
}



