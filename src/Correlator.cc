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





static int count_the_correlators = 1; 

static int count_the_zoomed_correlators = 1; 

static const UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance(); 



#ifndef RAD2DEG
#define RAD2DEG (180/M_PI)
#endif

#define PHI_SECTOR_ANGLE (360. / NUM_PHI)

//this should probably be defined elsewhere
#define ADU5_FORE_PHI 22.5



UCorrelator::Correlator::Correlator(int nphi, double phi_min, double phi_max, int ntheta, double theta_min, double theta_max)
  : hist(TString::Format("ucorr_corr_%d",count_the_correlators),"Correlator", nphi, phi_min, phi_max, ntheta, theta_min, theta_max), 
    norm(TString::Format("ucorr_norm_%d",count_the_correlators++),"Normalization", nphi, phi_min, phi_max, ntheta, theta_min, theta_max)
{
  
  
  trigcache = new TrigCache(nphi, (phi_max-phi_min)/nphi, phi_min, ntheta, (theta_max - theta_min)/ntheta,theta_min, ap); 

  disallowed_antennas = 0; 
  pad_factor = 1; 
  max_phi = 75; 

  memset(padded_waveforms, 0, NANTENNAS * sizeof(AnalysisWaveform*)); 
  memset(correlations, 0, NANTENNAS * NANTENNAS * sizeof(AnalysisWaveform*)); 
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

    for (int j = i; j < NANTENNAS; j++) 
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
  if (!padded_waveforms[ant1])
  {
//    printf("Copying and padding %d\n",ant1); 
     padded_waveforms[ant1] = new AnalysisWaveform(*ev->getFilteredGraph(ant1, pol)); 
     rms[ant1] = padded_waveforms[ant1]->even()->GetRMS(2); 
     padded_waveforms[ant1]->padEven(1); 
  }

  if (!padded_waveforms[ant2])
  {
//    printf("Copying and padding %d\n",ant2); 
      padded_waveforms[ant2] = new AnalysisWaveform(*ev->getFilteredGraph(ant2, pol)); 
//      printf("Computing rms!\n"); 
      rms[ant2] = padded_waveforms[ant2]->even()->GetRMS(2); 
      padded_waveforms[ant2]->padEven(1); 
  }


  if (!correlations[ant1][ant2])
  {
//    printf("Computing correlation %d %d\n", ant1, ant2); 
    correlations[ant1][ant2] = AnalysisWaveform::correlation(padded_waveforms[ant1],padded_waveforms[ant2],pad_factor, rms[ant1] * rms[ant2]); 
  }

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
  


  int closest[nant]; 


  ap->getClosestAntennas(phi, nant, closest, disallowed_antennas); 

  TrigCache cache(nphi, dphi, phi0, ntheta,dtheta,theta0, ap, true,nant, closest); 

  for (int ant_i = 0; ant_i < nant; ant_i++)
  {
    int ant1 = closest[ant_i]; 
    for (int ant_j = ant_i +1; ant_j < nant; ant_j++)
    {
      int ant2 = closest[ant_j]; 
      doAntennas(ant1, ant2, answer, &zoomed_norm, &cache); 
    }
  }

  return answer;
}


inline void UCorrelator::Correlator::doAntennas(int ant1, int ant2, TH2D * hist, 
                                                TH2I * norm, const TrigCache * cache = 0)
{
   int allowedFlag; 
   double lowerAngleThis, higherAngleThis, centerTheta1, centerTheta2, centerPhi1, centerPhi2;

   allowedFlag=allowedPhisPairOfAntennas(lowerAngleThis,higherAngleThis, 
                    centerTheta1, centerTheta2, centerPhi1, centerPhi2, 
                    ant1,ant2, max_phi, pol);

//   printf("HERE! allowedFlag:  %d\n", allowedFlag); 
   if(!allowedFlag) return; 

//   printf("Doing antennas (%d,%d)\n",ant1,ant2); 
//   printf("lowerAngleThis: %g higherAngleThis: %g\n", lowerAngleThis, higherAngleThis); 

   AnalysisWaveform * correlation = getCorrelation(ant1,ant2); 
   
   int nphibins = hist->GetNbinsX() + 2; 

   //find phi bin corresponding to lowerAngleThis and higherAngleThis

   int first_phi_bin = hist->GetXaxis()->FindFixBin(lowerAngleThis); 
   int last_phi_bin  = hist->GetXaxis()->FindFixBin(higherAngleThis); 

   bool must_wrap = (last_phi_bin < first_phi_bin) ; 

   for (int phibin = first_phi_bin; (phibin <= last_phi_bin) || must_wrap; phibin++)
   {
     if (must_wrap && phibin == nphibins-1)
     {
       phibin = 1; 
       must_wrap = false; 
     }


     double phi = cache? cache->phi[phibin-1] : hist->GetXaxis()->GetBinCenter(phibin);
     double dphi1 = FFTtools::wrap(phi - centerPhi1,360,0); 
     double dphi2 = FFTtools::wrap(phi - centerPhi2,360,0); 

     for (int thetabin = 1; thetabin <= hist->GetNbinsY(); thetabin++)
     {
       double theta = cache ? cache->theta[thetabin-1] : -hist->GetYaxis()->GetBinCenter(thetabin); //doh
       double dtheta1 = FFTtools::wrap(theta - centerTheta1,360,0); 
       double dtheta2 = FFTtools::wrap(theta - centerTheta2,360,0); 

       // check if in beam width 
       if (dphi1*dphi1 + dtheta1*dtheta1 > max_phi * max_phi) continue; 
       if (dphi2*dphi2 + dtheta2*dtheta2 > max_phi * max_phi) continue; 

       double timeExpected = cache? getDeltaTFast(ant1, ant2, phibin-1, thetabin-1,pol,cache)
                                  : getDeltaT(ant1, ant2, phi, theta,pol); 

       //add in off-axis antenna delay here
       if (groupDelayFlag)
       {
            Double_t delay1=getAntennaGroupDelay(dphi1,theta);
            Double_t delay2=getAntennaGroupDelay(dphi2,theta);
            timeExpected +=(delay1-delay2);
       }
         
       //TODO: add additional interpolation methods
       double val = correlation->evalEven(timeExpected); 

       int bin = phibin + thetabin * nphibins ; 


       hist->GetArray()[bin]+= val; 
       norm->GetArray()[bin]++;

     }
   }
}

void UCorrelator::Correlator::compute(const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t whichpol) 
{

  TStopwatch sw; 

  pol = whichpol; 
  ev = event; 
  hist.Reset(); 
  norm.Reset(); 

  reset(); 
  

  for (int ant1 = 0; ant1 < NANTENNAS; ant1++)
  {
    if (disallowed_antennas & (1 << ant1)) continue; 

    for (int ant2 = ant1+1; ant2 < NANTENNAS; ant2++) 
    {
      if (disallowed_antennas & (1 << ant2)) continue; 

      doAntennas(ant1, ant2, &hist, &norm, trigcache); 
    }
  }


  int nonzero = 0;
  //only keep values with at least 3 contributing antennas 
  for (int i = 0; i < (hist.GetNbinsX()+2) * (hist.GetNbinsY()+2); i++) 
  {
    double val = hist.GetArray()[i]; 
    if (val == 0) continue;
    int this_norm = norm.GetArray()[i]; 
    hist.GetArray()[i] = this_norm > 2 ? val/this_norm : 0;
  //  printf("%d %g %d\n",  i,  val, this_norm); 
    nonzero++; 
  }

  hist.SetEntries(nonzero); 
  norm.SetEntries(nonzero); 
  sw.Print("u");
}



UCorrelator::Correlator::~Correlator()
{

  delete trigcache; 
  reset(); 
}


void UCorrelator::Correlator::dumpDeltaTs(const char * fname) const
{

  TFile f(fname,"RECREATE"); 


  TTree * tree = new TTree("delays","Delays"); 

  int ant1, ant2; 
  double phi, theta, delta_t, group_delay; 

  tree->Branch("ant1",&ant1); 
  tree->Branch("ant2",&ant2); 
  tree->Branch("phi",&phi); 
  tree->Branch("theta",&theta); 
  tree->Branch("delta_t",&delta_t); 
  tree->Branch("group_delay",&group_delay); 


 

  for (ant1= 0; ant1 < NANTENNAS; ant1++) 
  {
    for (ant2 = ant1+1; ant2 < NANTENNAS; ant2++)
    {
      for (phi = -180; phi <=180; phi += 2)
      {
        for (theta = -90; theta <=90; theta +=2) 
        {
          delta_t = getDeltaT(ant1,ant2,phi,theta,pol); 

           if (groupDelayFlag)
           {
              double dphi1 = FFTtools::wrap(phi - ap->phiAnt[pol][ant1],360,0); 
              double dphi2 = FFTtools::wrap(phi - ap->phiAnt[pol][ant2],360,0); 
              double delay1=getAntennaGroupDelay(dphi1,theta);
              double delay2=getAntennaGroupDelay(dphi2,theta);
              group_delay = delay1-delay2; 
              delta_t +=group_delay;
           }
     
          f.cd(); 
          tree->Fill(); 
        }
      }
    }
  }

  f.cd(); 
  tree->Write(); 
  f.Close(); 


}
