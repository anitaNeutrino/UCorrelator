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
#endif

namespace UCorrelator
{
  class CorrelatorLocks
  {
    public:

#ifdef UCORRELATOR_OPENMP
      omp_lock_t waveform_locks[NANTENNAS]; 
      omp_lock_t correlation_locks[NANTENNAS][NANTENNAS]; 

      omp_lock_t hist_lock; 
      omp_lock_t norm_lock; 


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

        omp_init_lock(&hist_lock); 
        omp_init_lock(&norm_lock); 


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

        omp_destroy_lock(&hist_lock); 
        omp_destroy_lock(&norm_lock); 
      }
#endif
  };
}




static int count_the_correlators = 1; 

static int count_the_zoomed_correlators = 1; 




#ifndef RAD2DEG
#define RAD2DEG (180/M_PI)
#endif

#define PHI_SECTOR_ANGLE (360. / NUM_PHI)



UCorrelator::Correlator::Correlator(int nphi, double phi_min, double phi_max, int ntheta, double theta_min, double theta_max, bool use_center, bool scale_by_cos_theta, double baseline_weight)
  : scale_cos_theta(scale_by_cos_theta) , baselineWeight(baseline_weight)
{
  TString histname = TString::Format("ucorr_corr_%d",count_the_correlators);
  TString normname = TString::Format("ucorr_norm_%d",count_the_correlators++);
  hist = new TH2D(histname.Data(),"Correlator", nphi, phi_min, phi_max, ntheta, theta_min, theta_max); 
  norm = new TH2I(normname.Data(),"Normalization", nphi, phi_min, phi_max, ntheta, theta_min, theta_max);
    
  hist->GetXaxis()->SetTitle("#phi"); 
  hist->GetYaxis()->SetTitle("-#theta"); 
  
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
  TH2I zoomed_norm(TString::Format("zoomed_norm_%d",count_the_zoomed_correlators), "Zoomed Correlation Normalization", 
                    nphi, phi0,phi1, 
                    ntheta, theta0, theta1); 
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
    if (!nant && disallowed_antennas & (1 << ant1)) continue; 
    for (int ant_j = ant_i +1; ant_j < n2loop; ant_j++)
    {
      int ant2 = nant ? closest[ant_j] : ant_j; 
      if (!nant && disallowed_antennas & (1 << ant2)) continue; 

      pairs.push_back(std::pair<int,int>(ant1,ant2));
    }
  }
 
  unsigned nit = pairs.size(); 

#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for
#endif
  for (unsigned it = 0; it < nit; it++)
  {
     doAntennas(pairs[it].first, pairs[it].second, answer, &zoomed_norm, &cache, center_point); 
  }


  int nonzero = 0;
  //only keep values with at least  contributing antennas 
  for (int i = 0; i < (answer->GetNbinsX()+2) * (answer->GetNbinsY()+2); i++) 
  {
    double val = answer->GetArray()[i]; 
    if (val == 0) continue;
    int this_norm = zoomed_norm.GetArray()[i]; 
    answer->GetArray()[i] = this_norm > 0  ? val/this_norm : 0;
    nonzero++; 
  }

  answer->SetEntries(nonzero); 
 
  return answer;
}


static inline bool between(double phi, double low, double high)
{

  double diff_highlow = fmod(high - low + 360, 360); 
  double diff_philow = fmod(phi - low + 360, 360); 

  return diff_philow < diff_highlow; 
}


inline void UCorrelator::Correlator::doAntennas(int ant1, int ant2, TH2D * hist, 
                                                TH2I * norm, const TrigCache * cache , 
                                                const double * center_point )
{
   int allowedFlag; 
   double lowerAngleThis, higherAngleThis, centerTheta1, centerTheta2, centerPhi1, centerPhi2;

   allowedFlag=allowedPhisPairOfAntennas(lowerAngleThis,higherAngleThis, 
                    centerTheta1, centerTheta2, centerPhi1, centerPhi2, 
                    ant1,ant2, max_phi, pol);


   if(!allowedFlag) return; 

//   printf("lowerAngleThis: %g higherAngleThis: %g\n", lowerAngleThis, higherAngleThis); 
   // More stringent check if we have a center point
   if (center_point && !between(center_point[0], lowerAngleThis, higherAngleThis))  return; 

   AnalysisWaveform * correlation = getCorrelation(ant1,ant2); 
   
   int nphibins = hist->GetNbinsX() + 2; 

   //find phi bin corresponding to lowerAngleThis and higherAngleThis

   int first_phi_bin = center_point ? 1 : hist->GetXaxis()->FindFixBin(lowerAngleThis); 
   int last_phi_bin  = center_point ? hist->GetNbinsX() : hist->GetXaxis()->FindFixBin(higherAngleThis); 

   if (first_phi_bin == 0) first_phi_bin = 1; 
   if (last_phi_bin == hist->GetNbinsX()+1) last_phi_bin = hist->GetNbinsX(); 
   bool must_wrap = (last_phi_bin < first_phi_bin) ; 


   //So the maximum number of bins is going to be the total number of bins in the histogram. We probably won't fill all of them, 
   //but memory is cheap and std::vector is slow  

   int maxsize = hist->GetNbinsY() * hist->GetNbinsX(); 


   //This is bikeshedding, but allocate it all contiguosly 
   int * alloc = new int[3*maxsize]; 
   int * phibins = alloc;
   int * thetabins = alloc + maxsize; 
   int * bins_to_fill = alloc + 2 *maxsize; 

   int nbins_used =0; 


  
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

     int ny = hist->GetNbinsY(); 


     for (int thetabin = 1; thetabin <= ny; thetabin++)
     {
       double theta =  cache->theta[thetabin-1]; 
//       double theta4width = center_point ? center_point[1] : theta; 
       double dtheta1 = center_point ? 0 : FFTtools::wrap(theta- centerTheta1,360,0); 
       double dtheta2 = center_point ? 0 : FFTtools::wrap(theta- centerTheta2,360,0); 

       // check if in beam width 
       if (!center_point && dphi1*dphi1 + dtheta1*dtheta1 > max_phi * max_phi) continue; 
       if (!center_point && dphi2*dphi2 + dtheta2*dtheta2 > max_phi * max_phi) continue; 

       phibins[nbins_used] = phibin; 
       thetabins[nbins_used] = thetabin; 
       bins_to_fill[nbins_used] = phibin + thetabin * nphibins;
       nbins_used++; 
     }
   }


   double * dalloc = new double[2 *nbins_used]; 
   double * vals_to_fill = dalloc; 
   double * times_to_fill = dalloc + nbins_used; 

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

#ifdef UCORRELATOR_OPENMP
   omp_set_lock(&locks->hist_lock);
#endif
   for (int bi = 0; bi < nbins_used; bi++)
   {
       double val = vals_to_fill[bi]; 
       int bin = bins_to_fill[bi]; 
       hist->GetArray()[bin]+= val; 
   }
#ifdef UCORRELATOR_OPENMP
   omp_unset_lock(&locks->hist_lock); 
   omp_set_lock(&locks->norm_lock);
#endif   
   for (int bi = 0; bi < nbins_used; bi++)
   {
       int bin = bins_to_fill[bi]; 
       norm->GetArray()[bin]++;
   }
#ifdef UCORRELATOR_OPENMP
   omp_unset_lock(&locks->norm_lock); 
#endif
   delete [] alloc; 
   delete [] dalloc; 
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
    if (disallowed_antennas & (1 << ant1)) continue; 

    for (int ant2 = ant1+1; ant2 < NANTENNAS; ant2++)
    {
      if (disallowed_antennas & (1 << ant2)) continue; 

      pairs.push_back(std::pair<int,int>(ant1,ant2));;
    }

  }

  unsigned nit = pairs.size(); 

#ifdef UCORRELATOR_OPENMP
  #pragma omp parallel for 
#endif
  for (unsigned it = 0; it < nit; it++)
  {
     doAntennas(pairs[it].first, pairs[it].second, hist, norm, trigcache[v]); 
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
#endif 

}


void UCorrelator::Correlator::dumpDeltaTs(const char * fname) const
{

  TFile f(fname,"RECREATE"); 

  TTree * positions = new TTree("positions","Positions"); 


  TTree * tree = new TTree("delays","Delays"); 

  int ant1, ant2, pol; 
  double phi, theta, delta_t, group_delay; 

  tree->Branch("pol",&pol); 
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
  positions->Branch("pol",&pol); 

  const AntennaPositions * ap = AntennaPositions::instance(); 

  for (pol = 0; pol < 2; pol++)
  {
    for (ant1= 0; ant1 < NANTENNAS; ant1++) 
    {
      ant_phi = ap->phiAnt[pol][ant1]; 
      ant_r = ap->rAnt[pol][ant1]; 
      ant_z = ap->zAnt[pol][ant1]; 
      positions->Fill(); 

      for (ant2 = ant1+1; ant2 < NANTENNAS; ant2++)
      {
        for (phi = 0; phi <=360; phi += 2)
        {
          for (theta = -90; theta <=90; theta +=2) 
          {
             delta_t = getDeltaT(ant1,ant2,phi,theta,(AnitaPol::AnitaPol_t) pol, groupDelayFlag); 

             if (groupDelayFlag)
             {
                double dphi1 = FFTtools::wrap(phi - ap->phiAnt[pol][ant1],360,0); 
                double dphi2 = FFTtools::wrap(phi - ap->phiAnt[pol][ant2],360,0); 
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


