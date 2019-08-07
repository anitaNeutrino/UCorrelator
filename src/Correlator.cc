#include "FilteredAnitaEvent.h" 
#include "TString.h"
#include "AntennaPositions.h"
#include "DeltaT.h"
#include "TTree.h"
#include "TFile.h"
#include "TrigCache.h"
#include "FFTtools.h"
#include <assert.h>
#include "AnitaGeomTool.h"
#include "Correlator.h"
#include "AnalysisWaveform.h"
#include "TStopwatch.h" 

#ifdef UCORRELATOR_OPENMP
#include <omp.h>

#define SECTIONS _Pragma("omp parallel sections")
#define SECTION _Pragma("omp section") 

#else 

#define SECTIONS if(true) 
#define SECTION if(true) 


#endif

namespace UCorrelator
{
  class CorrelatorLocks
  {
    public:

#ifdef UCORRELATOR_OPENMP
      omp_lock_t waveform_locks[NANTENNAS]; 
      omp_lock_t correlation_locks[NANTENNAS][NANTENNAS]; 



      CorrelatorLocks() 
      {
        for (int i = 0; i < NANTENNAS; i++)
        {
          omp_init_lock(&waveform_locks[i]);

          for (int j = 0; j < NANTENNAS; j++) 
          {
            omp_init_lock(&correlation_locks[i][j]);
          }
        }



      }
      ~CorrelatorLocks() 
      {
        for (int i = 0; i < NANTENNAS; i++)
        {
          omp_destroy_lock(&waveform_locks[i]);

          for (int j = 0; j < NANTENNAS; j++) 
          {
            omp_destroy_lock(&correlation_locks[i][j]);
          }
        }

      }
#endif
  };
}

static int nthreads() 
{
#ifdef UCORRELATOR_OPENMP
  return omp_get_max_threads(); 
#else
  return 1;
#endif
}




static int gettid() 
{
#ifdef UCORRELATOR_OPENMP
  return omp_get_thread_num(); 
#else
  return 0;
#endif
}




static int count_the_correlators = 1; 

static int count_the_zoomed_correlators = 1; 



/*Add all histograms to the first one ... quickly... assuming they're the same size and type */
template<class T,class A> 
static void combineHists(int N, T ** hists )
{
  if (N <2) return; 
  int nbins = (hists[0]->GetNbinsX() +2) * (hists[0]->GetNbinsY()+2); 
  __restrict A* sum = hists[0]->GetArray(); 
  for (int i = 1; i <N; i++ )
  {
    __restrict A* v = hists[i]->GetArray(); 
#pragma omp simd 
    for(int j = 0; j < nbins; j++) 
    {
      sum[j]+=v[j]; 
    }
  }
}



#ifndef RAD2DEG
#define RAD2DEG (180/M_PI)
#endif

#define PHI_SECTOR_ANGLE (360. / NUM_PHI)



UCorrelator::Correlator::Correlator(int nphi, double phi_min, double phi_max, int ntheta, double theta_min, double theta_max, bool use_center, bool scale_by_cos_theta, double baseline_weight, double gain_sigma)
  : scale_cos_theta(scale_by_cos_theta) , baselineWeight(baseline_weight), gainSigma(gain_sigma) 
{
  TString histname = TString::Format("ucorr_corr_%d",count_the_correlators);
  TString normname = TString::Format("ucorr_norm_%d",count_the_correlators++);

  hist = new TH2D(histname.Data(),"Correlator", nphi, phi_min, phi_max, ntheta, theta_min, theta_max); 
  hist->SetDirectory(0); 
  norm = new TH2D(normname.Data(),"Normalization", nphi, phi_min, phi_max, ntheta, theta_min, theta_max);
  norm->SetDirectory(0); 

  hist->GetXaxis()->SetTitle("#phi"); 
  hist->GetYaxis()->SetTitle("-#theta"); 

  hists.resize(nthreads()); 
  norms.resize(nthreads()); 
  hists[0] = hist;
  norms[0] = norm;

//  printf(":: %p %p\n",hists[0],norms[0]); 
  for (int i= 1; i < nthreads();i++)
  {
    hists[i] = (TH2D*) hist->Clone(TString(hist->GetName()) + TString::Format("_%d",i)); 
    hists[i]->SetDirectory(0); 
    norms[i] = (TH2D*) norm->Clone(TString(norm->GetName()) + TString::Format("_%d",i)); 
    norms[i]->SetDirectory(0); 
//    printf(":: %p %p\n",hists[i],norms[i]); 
  }



  use_bin_center = use_center;
  memset(trigcache,0, sizeof(trigcache)); 


  disallowed_antennas = 0; 
  pad_factor = 3; 
  max_phi = 75; 
  max_phi2 = max_phi*max_phi;

  memset(padded_waveforms, 0, NANTENNAS * sizeof(AnalysisWaveform*)); 
  memset(correlations, 0, NANTENNAS * NANTENNAS * sizeof(AnalysisWaveform*)); 

#ifdef UCORRELATOR_OPENMP
  locks = new CorrelatorLocks; 
#endif
  ev = 0; 
  groupDelayFlag = 1; 
}

static int allowedPhisPairOfAntennas(double &lowerAngle, double &higherAngle, double &centerTheta1, double &centerTheta2, double &centerPhi1, double &centerPhi2, int ant1, int ant2, double max_phi, AnitaPol::AnitaPol_t pol)
{

  int phi1=AnitaGeomTool::Instance()->getPhiFromAnt(ant1);
  int phi2=AnitaGeomTool::Instance()->getPhiFromAnt(ant2);
  int allowedFlag=0;
  
  int upperlimit=phi2+2;//2 phi sectors on either side
  int lowerlimit=phi2-2;

  if(upperlimit>NUM_PHI-1)upperlimit-=NUM_PHI;
  if(lowerlimit<0)lowerlimit+=NUM_PHI;

  if (upperlimit>lowerlimit){
    if (phi1<=upperlimit && phi1>=lowerlimit){//within 2 phi sectors of eachother
      allowedFlag=1;
    }
  }
  if (upperlimit<lowerlimit){
    if (phi1<=upperlimit || phi1>=lowerlimit){
      allowedFlag=1;

    }
  }
  
  double centerAngle1, centerAngle2;

  const UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance(); 

  if (allowedFlag==1)
  {
    centerAngle1=ap->phiAnt[pol][ant1]; 
    centerAngle2=ap->phiAnt[pol][ant2]; 
//    assert(centerAngle1 == ap->phiAnt[0][ant1]); 

    if (centerAngle2>centerAngle1)
    {
      lowerAngle=centerAngle2-max_phi;
      higherAngle=centerAngle1+max_phi;
    }
    else
    {
      lowerAngle=centerAngle1-max_phi;
      higherAngle=centerAngle2+max_phi; 
    }

    if (lowerAngle<0) lowerAngle+=360;
    if (higherAngle>360) higherAngle-=360;
    
  }
  else
 {

    centerAngle1= 0; 
    centerAngle2= 0; 
  }
  centerTheta1=10;//degrees down
  centerTheta2=10;//degrees down
  centerPhi1=centerAngle1;
  centerPhi2=centerAngle2;
  
  return allowedFlag;

}

void UCorrelator::Correlator::reset() 
{
for (int i = 0; i < NANTENNAS; i++)
  {
    if (padded_waveforms[i]) 
    {
      delete padded_waveforms[i]; 
      padded_waveforms[i] = 0; 
    }

    for (int j = 0; j < NANTENNAS; j++) 
    {
      if (correlations[i][j])
      {
        delete correlations[i][j]; 
        correlations[i][j] = 0; 
      }
    }
  }


}




AnalysisWaveform * UCorrelator::Correlator::getCorrelation(int ant1, int ant2) 
{
//  printf("%d %d / %d \n",ant1,ant2, NANTENNAS); 


#ifdef UCORRELATOR_OPENMP
#ifndef FFTTOOLS_COMPILED_WITH_OPENMP
#pragma omp critical (getCorrelation)
  {
#endif
#endif


#ifdef UCORRELATOR_OPENMP
  omp_set_lock(&locks->waveform_locks[ant1]); 
#endif
  if (!padded_waveforms[ant1])
  {
//    printf("Copying and padding %d\n",ant1); 
     padded_waveforms[ant1] = new AnalysisWaveform(*ev->getFilteredGraph(ant1, pol)); 
     rms[ant1] = padded_waveforms[ant1]->even()->GetRMS(2); 
     padded_waveforms[ant1]->padEven(1); 
  }
#ifdef UCORRELATOR_OPENMP
  omp_unset_lock(&locks->waveform_locks[ant1]); 
  omp_set_lock(&locks->waveform_locks[ant2]); 
#endif

  if (!padded_waveforms[ant2])
  {
//    printf("Copying and padding %d\n",ant2); 
      padded_waveforms[ant2] = new AnalysisWaveform(*ev->getFilteredGraph(ant2, pol)); 
//      printf("Computing rms!\n"); 
      rms[ant2] = padded_waveforms[ant2]->even()->GetRMS(2); 
      padded_waveforms[ant2]->padEven(1); 
  }

#ifdef UCORRELATOR_OPENMP
  omp_unset_lock(&locks->waveform_locks[ant2]); 
  omp_set_lock(&locks->correlation_locks[ant1][ant2]); 
#endif

  if (!correlations[ant1][ant2])
  {
#ifdef UCORRELATOR_OPENMP
    omp_set_lock(&locks->waveform_locks[ant1]); 
    omp_set_lock(&locks->waveform_locks[ant2]); 
#endif
//    printf("Computing correlation %d %d\n", ant1, ant2); 
    correlations[ant1][ant2] = AnalysisWaveform::correlation(padded_waveforms[ant1],padded_waveforms[ant2],pad_factor, rms[ant1] * rms[ant2]); 

#ifdef UCORRELATOR_OPENMP
    omp_unset_lock(&locks->waveform_locks[ant2]); 
    omp_unset_lock(&locks->waveform_locks[ant1]); 
#endif
  }


#ifdef UCORRELATOR_OPENMP
  omp_unset_lock(&locks->correlation_locks[ant1][ant2]); 
#ifndef FFTTOOLS_COMPILED_WITH_OPENMP
  }
#endif
#endif
  return correlations[ant1][ant2]; 
}



TH2D * UCorrelator::Correlator::computeZoomed(double phi, double theta, int nphi, double dphi, int ntheta, double dtheta, int nant, TH2D * answer) 
{

  if (!ev) 
  {
    fprintf(stderr, "Must call Correlator::compute() prior to Correlator::computeZoomed!!!!"); 
    return 0; 
  }

  double phi0 = phi - dphi * nphi/2; 
  double phi1 = phi + dphi * nphi/2; 
  double theta0 = theta - dtheta * ntheta/2; 
  double theta1 = theta + dtheta * ntheta/2; 
  TH2D* zoomed_norm = new TH2D(TString::Format("zoomed_norm_%d",count_the_zoomed_correlators), "Zoomed Correlation Normalization", 
                    nphi, phi0,phi1, 
                    ntheta, theta0, theta1); 
  zoomed_norm->SetDirectory(0); 

  if (answer) 
  {
    answer->SetBins(nphi, phi0, phi1, 
                    ntheta, theta0, theta1); 
    answer->Reset(); 
  }
  else
  {
    answer = new TH2D(TString::Format("zoomed_corr_%d", count_the_zoomed_correlators++), "Zoomed Correlation", 
                    nphi, phi0, phi1, 
                    ntheta, theta0, theta1); 
    answer->SetDirectory(0); 
  }
  
  answer->GetXaxis()->SetTitle("#phi"); 
  answer->GetYaxis()->SetTitle("-#theta"); 

  const AntennaPositions * ap = AntennaPositions::instance(); 

  int closest[nant]; // is it a problem if nant is 0? 
  if (nant) 
  {
    memset(closest,0,sizeof(closest)); 
    nant = ap->getClosestAntennas(phi, nant, closest, disallowed_antennas); 
  }

  TrigCache cache(nphi, dphi, phi0, ntheta,dtheta,theta0, ap, true,nant, nant ? closest : 0); 

  int n2loop = nant ? nant : NANTENNAS;  

  std::vector<std::pair<int,int> > pairs; 
  pairs.reserve(n2loop); 

  double center_point[2];
  center_point[0] = phi; 
  center_point[1] = theta; 

  for (int ant_i = 0; ant_i < n2loop; ant_i++)
  {
    int ant1 = nant ? closest[ant_i] : ant_i; 
    if (!nant && disallowed_antennas & (1ul << ant1)) continue; 
    for (int ant_j = ant_i +1; ant_j < n2loop; ant_j++)
    {
      int ant2 = nant ? closest[ant_j] : ant_j; 
      if (!nant && disallowed_antennas & (1ul << ant2)) continue; 

      pairs.push_back(std::pair<int,int>(ant1,ant2));
    }
  }
 
  unsigned nit = pairs.size(); 

 

  /* lock contention for the hist / norm lock is killing parallelization */ 
  std::vector<TH2D*> zoomed_hists(nthreads()); 
  std::vector<TH2D*> zoomed_norms(nthreads()); 

  zoomed_hists[0] = answer;
  zoomed_norms[0] = zoomed_norm;

//  printf(":: %p %p\n",hists[0],norms[0]); 
  for (int i= 1; i < nthreads();i++)
  {
    zoomed_hists[i] = (TH2D*) answer->Clone(TString(answer->GetName()) + TString::Format("_%d",i)); 
    zoomed_hists[i]->SetDirectory(0); 
    zoomed_norms[i] = (TH2D*) zoomed_norm->Clone(TString(norm->GetName()) + TString::Format("_%d",i)); 
    norms[i]->SetDirectory(0); 
//  zoomed_  printf(":: %p %p\n",hists[i],norms[i]); 
  }


#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for
#endif
  for (unsigned it = 0; it < nit; it++)
  {
     doAntennas(pairs[it].first, pairs[it].second, &zoomed_hists[0], &zoomed_norms[0], &cache, center_point); 
  }


  /* combine histograms */ 
  
SECTIONS
  {
SECTION
    combineHists<TH2D,double>(nthreads(), &zoomed_hists[0]); 

SECTION
    combineHists<TH2D,double>(nthreads(), &zoomed_norms[0]); 
  }

  for (int i =1; i < nthreads(); i++) 
  {
    delete zoomed_hists[i]; 
    delete zoomed_norms[i]; 
  }


  int nonzero = 0;
  //only keep values with contributing antennas 
  for (int i = 0; i < (answer->GetNbinsX()+2) * (answer->GetNbinsY()+2); i++) 
  {
    double val = answer->GetArray()[i]; 
    if (val == 0) continue;
    double this_norm = zoomed_norm->GetArray()[i]; 
    answer->GetArray()[i] = this_norm > 0  ? val/this_norm : 0;
    nonzero++; 
  }

  answer->SetEntries(nonzero); 
 
  delete zoomed_norm; 
  return answer;
}


static inline bool between(double phi, double low, double high)
{

  double diff_highlow = fmod(high - low + 360, 360); 
  double diff_philow = fmod(phi - low + 360, 360); 

  return diff_philow < diff_highlow; 
}


inline void UCorrelator::Correlator::doAntennas(int ant1, int ant2, TH2D ** these_hists, 
                                                TH2D ** these_norms, const UCorrelator::TrigCache * cache , 
                                                const double * center_point )
{
   int allowedFlag; 
   double lowerAngleThis, higherAngleThis, centerTheta1, centerTheta2, centerPhi1, centerPhi2;

   allowedFlag=allowedPhisPairOfAntennas(lowerAngleThis,higherAngleThis, 
                    centerTheta1, centerTheta2, centerPhi1, centerPhi2, 
                    ant1,ant2, max_phi, pol);


   assert(ant1 < 48); 
   assert(ant2 < 48); 
   if(!allowedFlag) return; 

   TH2D * the_hist  = these_hists[gettid()]; 
   TH2D * the_norm  = these_norms[gettid()]; 

//   printf("lowerAngleThis: %g higherAngleThis: %g\n", lowerAngleThis, higherAngleThis); 
   // More stringent check if we have a center point
   if (center_point && !between(center_point[0], lowerAngleThis, higherAngleThis))  return; 

   AnalysisWaveform * correlation = getCorrelation(ant1,ant2); 
   
   int nphibins = the_hist->GetNbinsX() + 2; 

   //find phi bin corresponding to lowerAngleThis and higherAngleThis

   int first_phi_bin = center_point ? 1 : the_hist->GetXaxis()->FindFixBin(lowerAngleThis); 
   int last_phi_bin  = center_point ? the_hist->GetNbinsX() : the_hist->GetXaxis()->FindFixBin(higherAngleThis); 

   if (first_phi_bin == 0) first_phi_bin = 1; 
   if (last_phi_bin == the_hist->GetNbinsX()+1) last_phi_bin = the_hist->GetNbinsX(); 
   bool must_wrap = (last_phi_bin < first_phi_bin) ; 


   //So the maximum number of bins is going to be the total number of bins in the histogram. We probably won't fill all of them, 
   //but memory is cheap and std::vector is slow  

   int maxsize = the_hist->GetNbinsY() * the_hist->GetNbinsX(); 

   double phi_diff1 = FFTtools::wrap(centerPhi1 - centerPhi2, 360, 0); 
   double phi_diff2 = FFTtools::wrap(centerPhi2 - centerPhi1, 360, 0); 
   double baseline_phi = (fabs(phi_diff1) < fabs(phi_diff2)) ? centerPhi2 + phi_diff1 / 2 : centerPhi1 + phi_diff2 / 2; 
   baseline_phi = FFTtools::wrap(baseline_phi, 360); 

   double baseline_theta = (centerTheta1 + centerTheta2) / 2;
 

   //This is bikeshedding, but allocate it all contiguosly 
   int * alloc = new int[3*maxsize]; 
   int * phibins = alloc;
   int * thetabins = alloc + maxsize; 
   int * bins_to_fill = alloc + 2 *maxsize; 

   int nbins_used =0; 

   double * gain_weights = 0;
   if (gainSigma && !center_point) gain_weights = new double[maxsize]; 

  
   for (int phibin = first_phi_bin; (phibin <= last_phi_bin) || must_wrap; phibin++)
   {
     if (must_wrap && phibin == nphibins-1)
     {
       phibin = 1; 
       must_wrap = false; 
     }
     
     double phi =  cache->phi[phibin-1] ; 
//     double phi4width = center_point ? center_point[0] : phi; 
     double dphi1 = center_point ? 0 : FFTtools::wrap(phi - centerPhi1,360,0); 
     double dphi2 = center_point ? 0 : FFTtools::wrap(phi - centerPhi2,360,0); 


     //Check if in beam width in phi 
     if (!center_point && fabs(dphi1)  > max_phi) continue; 
     if (!center_point && fabs(dphi2)  > max_phi) continue; 

     int ny = the_hist->GetNbinsY(); 


     for (int thetabin = 1; thetabin <= ny; thetabin++)
     {
       double theta =  cache->theta[thetabin-1]; 
//       double theta4width = center_point ? center_point[1] : theta; 
       double dtheta1 = center_point ? 0 : FFTtools::wrap(theta - centerTheta1, 360, 0); 
       double dtheta2 = center_point ? 0 : FFTtools::wrap(theta - centerTheta2, 360, 0); 

       // check if in beam width 
       if (!center_point && dphi1*dphi1 + dtheta1*dtheta1 > max_phi * max_phi) continue; 
       if (!center_point && dphi2*dphi2 + dtheta2*dtheta2 > max_phi * max_phi) continue; 

       phibins[nbins_used] = phibin; 
       thetabins[nbins_used] = thetabin; 
       bins_to_fill[nbins_used] = phibin + thetabin * nphibins;
       if (gainSigma && !center_point) 
       {
         //Matt Mottram weighted by the baseline angle difference somehow
         
        double Dphi = FFTtools::wrap(phi - baseline_phi, 360, 0); 
        double Dtheta = theta - baseline_theta;

        gain_weights[nbins_used] = TMath::Gaus(TMath::Sqrt(Dphi * Dphi + Dtheta * Dtheta), 0, gainSigma);
       }
       nbins_used++; 
     }
   }


   double * dalloc = new double[2 *nbins_used]; 
   double * __restrict__ vals_to_fill = dalloc; 
   double * __restrict__ times_to_fill = dalloc + nbins_used; 

  //TODO vectorize this
   for (int i = 0; i < nbins_used; i++)
   {
       
       int phibin = phibins[i];; 
       int thetabin = thetabins[i]; 
       times_to_fill[i] = getDeltaTFast(ant1, ant2, phibin-1, thetabin-1,pol,cache, groupDelayFlag); 
   }

   correlation->evalEven(nbins_used, times_to_fill, vals_to_fill); 


   if (scale_cos_theta)
   {
     for(int i = 0; i < nbins_used; i++)
     {
       vals_to_fill[i] *= cache->cos_theta[thetabins[i]-1];
     }
   }


   if (baselineWeight)
   {
     double wgt = TMath::Power(cache->ap->distance(ant1, ant2, pol), baselineWeight); 
     for (int i = 0; i < nbins_used; i++) 
     {
       vals_to_fill[i] *= wgt; 
     }
   }



   double * __restrict__ the_arr = the_hist->GetArray();
   double * __restrict__ the_norm_arr = the_norm->GetArray();

   for (int bi = 0; bi < nbins_used; bi++)
   {
       double val = vals_to_fill[bi]; 
       int bin = bins_to_fill[bi]; 

       if (gainSigma && !center_point)
       {
         the_arr[bin] += val * gain_weights[bi];
         the_norm_arr[bin] += gain_weights[bi];
       }
       else
       {
         the_arr[bin]+=  val; 
         the_norm_arr[bin]+= 1;
       }
   }


   delete [] alloc; 
   delete [] dalloc; 
   if (gain_weights) delete [] gain_weights; 
}


void UCorrelator::Correlator::compute(const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t whichpol) 
{

//  TStopwatch sw; 

  pol = whichpol; 
  ev = event; 
  hist->Reset(); 
  norm->Reset(); 
  reset(); 

  //alright, we have to be able to dispatch between different versions 
  //because who knows what crazy things we might be asked to do. 
  // I suppose we could just require the user to make a different Correlator for each ANITA version... but whatever
  
  int v = event->getAnitaVersion(); 
  AnitaVersion::set(v); 

  if (! trigcache[v])
  {
    int nphi = hist->GetNbinsX(); 
    double phi_max = hist->GetXaxis()->GetXmax(); 
    double phi_min = hist->GetXaxis()->GetXmin(); 
    int ntheta = hist->GetNbinsY(); 
    double theta_max = hist->GetYaxis()->GetXmax(); 
    double theta_min = hist->GetYaxis()->GetXmin(); 
    trigcache[v] = new TrigCache(nphi, (phi_max-phi_min)/nphi, phi_min, ntheta, (theta_max - theta_min)/ntheta,theta_min, UCorrelator::AntennaPositions::instance(v), use_bin_center); 
  }


  //precompute antenna combinations 


  std::vector<std::pair<int,int> > pairs; 
  pairs.reserve(NANTENNAS *NANTENNAS/2); 

  for (int ant1 = 0; ant1 < NANTENNAS; ant1++)
  {
    if (disallowed_antennas & (1ul << ant1)) continue; 

    for (int ant2 = ant1+1; ant2 < NANTENNAS; ant2++)
    {
      if (disallowed_antennas & (1ul << ant2)) continue; 

      pairs.push_back(std::pair<int,int>(ant1,ant2));;
    }

  }

  unsigned nit = pairs.size(); 

  for (int i= 1; i < nthreads();i++)
  {
    hists[i]->Reset(); 
    norms[i]->Reset();
  }

#ifdef UCORRELATOR_OPENMP
  #pragma omp parallel for 
#endif
  for (unsigned it = 0; it < nit; it++)
  {
     doAntennas(pairs[it].first, pairs[it].second, &hists[0], &norms[0], trigcache[v]); 
  }

  /* combine histograms */ 
SECTIONS
{

  SECTION
    combineHists<TH2D,double>(nthreads(),&hists[0]); 
  SECTION
    combineHists<TH2D,double>(nthreads(),&norms[0]); 
}


  int nonzero = 0;
  //only keep values with at least 3 contributing antennas 
  for (int i = 0; i < (hist->GetNbinsX()+2) * (hist->GetNbinsY()+2); i++) 
  {
    double val = hist->GetArray()[i]; 
    if (val == 0) continue;
    int this_norm = norm->GetArray()[i]; 
    hist->GetArray()[i] = this_norm > 2 ? val/this_norm : 0;

  //  printf("%d %g %d\n",  i,  val, this_norm); 
    nonzero++; 
  }

  hist->SetEntries(nonzero); 
  norm->SetEntries(nonzero); 
//  sw.Print("u");
}



UCorrelator::Correlator::~Correlator()
{

  for (int i = 0; i < NUM_ANITAS+1; i++) 
  {
    if (trigcache[i])
    {
      delete trigcache[i]; 
    }
  }

  reset(); 
  delete hist; 
  delete norm; 

#ifdef UCORRELATOR_OPENMP
  delete locks; 
  for (int i = 1; i < nthreads(); i++)
  {
    delete hists[i];
    delete norms[i];

  }
#endif 

}


void UCorrelator::Correlator::dumpDeltaTs(const char * fname) const
{

  TFile f(fname,"RECREATE"); 

  TTree * positions = new TTree("positions","Positions"); 


  TTree * tree = new TTree("delays","Delays"); 

  int ant1, ant2, ipol; 
  double phi, theta, delta_t, group_delay; 

  tree->Branch("pol",&ipol); 
  tree->Branch("ant1",&ant1); 
  tree->Branch("ant2",&ant2); 
  tree->Branch("phi",&phi); 
  tree->Branch("theta",&theta); 
  tree->Branch("delta_t",&delta_t); 
  tree->Branch("group_delay",&group_delay); 

  double ant_phi, ant_r, ant_z; 
  positions->Branch("ant",&ant1); 
  positions->Branch("phi",&ant_phi); 
  positions->Branch("r",&ant_r); 
  positions->Branch("z",&ant_z); 
  positions->Branch("pol",&ipol); 

  const AntennaPositions * ap = AntennaPositions::instance(); 

  for (ipol = 0; ipol < 2; ipol++)
  {
    for (ant1= 0; ant1 < NANTENNAS; ant1++) 
    {
      ant_phi = ap->phiAnt[ipol][ant1]; 
      ant_r = ap->rAnt[ipol][ant1]; 
      ant_z = ap->zAnt[ipol][ant1]; 
      positions->Fill(); 

      for (ant2 = ant1+1; ant2 < NANTENNAS; ant2++)
      {
        for (phi = 0; phi <=360; phi += 2)
        {
          for (theta = -90; theta <=90; theta +=2) 
          {
             delta_t = getDeltaT(ant1,ant2,phi,theta,(AnitaPol::AnitaPol_t) ipol, groupDelayFlag); 

             if (groupDelayFlag)
             {
                double dphi1 = FFTtools::wrap(phi - ap->phiAnt[ipol][ant1],360,0); 
                double dphi2 = FFTtools::wrap(phi - ap->phiAnt[ipol][ant2],360,0); 
                double delay1=getAntennaGroupDelay(dphi1,theta);
                double delay2=getAntennaGroupDelay(dphi2,theta);
                group_delay = delay1-delay2; 
             }
             f.cd(); 
             tree->Fill(); 
          }
        }
      }
    }
  }

  f.cd(); 
  positions->Write(); 
  tree->Write(); 
  f.Close(); 


}


