#include "FilteredAnitaEvent.h" 
#include "TString.h"
#include "AntennaPositions.h"
#include "FFTtools.h"
#include "AnitaGeomTool.h"
#include "Correlator.h"
#include "AnalysisWaveform.h"


static int count_the_correlators = 1; 

static int count_the_zoomed_correlators = 1; 

static UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance(); 


#define DEG2RAD (M_PI/180)
#define RAD2DEG (180/M_PI)

#define PHI_SECTOR_ANGLE (360 / PHI_SECTORS)

//this should probably be dfined elsewher
#define ADU5_FORE_PHI 22.5



UCorrelator::Correlator::Correlator(int nphi, double phi_min, double phi_max, int ntheta, double theta_min, double theta_max)
  : hist(TString::Format("ucorr_corr_%d",count_the_correlators),"Correlator", nphi, phi_min, phi_max, ntheta, theta_min, theta_max), 
    norm(TString::Format("ucorr_norm_%d",count_the_correlators++),"Normalization", nphi, phi_min, phi_max, ntheta, theta_min, theta_max)
{
  
  
  phi_table = new double[nphi]; 
  theta_table = new double[ntheta]; 
  cos_phi_table = new double[nphi]; 
  tan_theta_table = new double[ntheta]; 
  cos_theta_table = new double[ntheta]; 
  compute_table(); 

  pad_factor = 1; 
  max_phi = 75; 

  memset(padded_waveforms, 0, NANTENNAS * sizeof(AnalysisWaveform*)); 
  memset(correlations, 0, NANTENNAS * NANTENNAS * sizeof(AnalysisWaveform*)); 
  ev = 0; 
  groupDelayFlag = 1; 
}

void UCorrelator::Correlator::compute_table()
{
  double phi_min = hist.GetXaxis()->GetXmin(); 
  double theta_min = hist.GetYaxis()->GetXmin(); 
  double phi_max = hist.GetXaxis()->GetXmax(); 
  double theta_max = hist.GetYaxis()->GetXmax(); 
  int nphi = hist.GetNbinsX(); 
  int ntheta = hist.GetNbinsY(); 

  double phi0 = phi_min * DEG2RAD;  
  double dphi = (phi_max-phi_min)/nphi * DEG2RAD;  

  double theta0 = theta_min * DEG2RAD;  
  double dtheta = (theta_max-theta_min)/ntheta * DEG2RAD;  

  for (int i = 0; i < nphi; i++) 
  {
    double phi = phi0 + (i+0.5) * dphi; 
    phi_table[i] = phi;
    cos_phi_table[i] = cos(phi); 
  }

  for (int i = 0; i < ntheta; i++) 
  {
    double theta = theta0 + (i+0.5) * dtheta; 
    theta_table[i] = theta; 
    cos_theta_table[i] = cos(theta); 
    tan_theta_table[i] = tan(theta); 
  }
}

static int allowedPhisPairOfAntennas(double &lowerAngle, double &higherAngle, double &centerTheta1, double &centerTheta2, double &centerPhi1, double &centerPhi2, int ant1, int ant2, double max_phi)
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
    centerAngle1=phi1*PHI_SECTOR_ANGLE-ADU5_FORE_PHI;
    centerAngle2=phi2*PHI_SECTOR_ANGLE-ADU5_FORE_PHI;

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
  if (!padded_waveforms[ant1])
  {
     padded_waveforms[ant1] = new AnalysisWaveform(*ev->getFilteredGraph(ant1, pol)); 
     rms[ant1] = padded_waveforms[ant1]->even()->GetRMS(2); 
     padded_waveforms[ant1]->padEven(1); 
  }

  if (!padded_waveforms[ant2])
  {
      padded_waveforms[ant2] = new AnalysisWaveform(*ev->getFilteredGraph(ant2, pol)); 
      rms[ant2] = padded_waveforms[ant2]->even()->GetRMS(2); 
      padded_waveforms[ant2]->padEven(1); 
  }


  if (!correlations[ant1][ant2])
  {
    correlations[ant1][ant2] = AnalysisWaveform::correlation(padded_waveforms[ant1],padded_waveforms[ant2],pad_factor, rms[ant1] * rms[ant2]); 
  }

  return correlations[ant1][ant2]; 
}



TH2D * UCorrelator::Correlator::computeZoomed(double phi, double theta, int ntheta, double dtheta, int nphi, double dphi, int nant, TH2D * answer) 
{

  if (!ev) 
  {
    fprintf(stderr, "Must call Correlator::compute() prior to Correlator::computeZoomed!!!!"); 
    return 0; 
  }

  TH2I norm(TString::Format("zoomed_norm_%d",count_the_zoomed_correlators), "Zoomed Correlation Normalization", 
                    nphi, phi-dphi*nphi/2, phi + dphi * nphi/2, 
                    ntheta, theta - dtheta *ntheta/2, theta + dtheta * ntheta/2); 
  if (answer) 
  {
    answer->SetBins(nphi, phi-dphi*nphi/2, phi + dphi * nphi/2, 
                    ntheta, theta - dtheta *ntheta/2, theta + dtheta * ntheta/2); 
    answer->Reset(); 
  }
  else
  {
    answer = new TH2D(TString::Format("zoomed_corr_%d", count_the_zoomed_correlators++), "Zoomed Correlation", 
                    nphi, phi-dphi*nphi/2, phi + dphi * nphi/2, 
                    ntheta, theta - dtheta *ntheta/2, theta + dtheta * ntheta/2); 
  }
  


  int closest[nant]; 


  ap->getClosestAntennas(phi, nant, closest, disallowed_antennas); 


  for (int ant_i = 0; ant_i < nant; ant_i++)
  {
    int ant1 = closest[ant_i]; 
    for (int ant_j = ant_i +1; ant_j < nant; ant_j++)
    {
      int ant2 = closest[ant_j]; 
      doAntennas(ant1, ant2, &hist, &norm, false); 
    }
  }

  return answer;
}


inline void UCorrelator::Correlator::doAntennas(int ant1, int ant2, TH2D * hist, 
                                                TH2I * norm, bool useArray)
{
   int allowedFlag; 
   double lowerAngleThis, higherAngleThis, centerTheta1, centerTheta2, centerPhi1, centerPhi2;
   allowedFlag=allowedPhisPairOfAntennas(lowerAngleThis,higherAngleThis, 
                    centerTheta1, centerTheta2, centerPhi1, centerPhi2, 
                    ant1,ant2, max_phi);

   if(!allowedFlag) return; 

   AnalysisWaveform * correlation = getCorrelation(ant1,ant2); 
   int nphibins = hist->GetNbinsX(); 

   //find phi bin corresponding to lowerAngleThis and higherAngleThis

   int first_phi_bin = hist->GetXaxis()->FindFixBin(lowerAngleThis); 
   int last_phi_bin = hist->GetXaxis()->FindFixBin(higherAngleThis); 

   for (int phibin = first_phi_bin; phibin <= last_phi_bin; phibin++)
   {
     double phi = useArray ? phi_table[phibin] : hist->GetXaxis()->GetBinCenter(phibin);
     double dphi1 = FFTtools::wrap(phi - centerPhi1,360,0); 
     double dphi2 = FFTtools::wrap(phi - centerPhi2,360,0); 

     for (int thetabin = 1; thetabin <= hist->GetNbinsY(); thetabin++)
     {
       double theta = useArray ? theta_table[thetabin] : hist->GetYaxis()->GetBinCenter(thetabin); 

       double dtheta1 = FFTtools::wrap(theta - centerTheta1,360,0); 
       double dtheta2 = FFTtools::wrap(theta - centerTheta2,360,0); 

       // check if in beam width 
       if (dphi1*dphi1 + dtheta1 *dtheta1 > max_phi * max_phi) continue; 
       if (dphi2*dphi2 + dtheta2 *dtheta2 > max_phi * max_phi) continue; 

       double timeExpected  = useArray ? getDeltaTFast(ant1, ant2, phibin, thetabin) : getDeltaT(ant1, ant2, phi, theta); 

       //add in off-axis antenna delay here
       if (groupDelayFlag)
       {
            Double_t delay1=getGroupDelay(dphi1,theta);
            Double_t delay2=getGroupDelay(dphi2,theta);
            timeExpected +=(delay1-delay2);
       }
         
       //TODO: add additional interpolation methods
       double val = correlation->even()->Eval(timeExpected); 

       int bin = phibin + thetabin  * (nphibins + 2); 

       hist->GetArray()[bin] += val; 
       norm->GetArray()[bin]++ ;

     }
   }
}

void UCorrelator::Correlator::compute(const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t pol) 
{
  ev = event; 

  hist.Reset(); 
  norm.Reset(); 


  for (int ant1 = 0; ant1 < NANTENNAS; ant1++)
  {
    if (disallowed_antennas & (1 << ant1)) continue; 

    for (int ant2 = ant1+1; ant2 < NANTENNAS; ant2++) 
    {
      if (disallowed_antennas & (1 << ant2)) continue; 

      doAntennas(ant1, ant2, &hist, &norm, true); 
    }
  }


  //only keep values with at least 3 contributing antennas 
  for (int i = 0; i < (hist.GetNbinsX()+2) * (hist.GetNbinsY()+2); i++) 
  {
    double val = hist.GetArray()[i]; 
    if (val == 0) continue;
    int this_norm = norm.GetArray()[i]; 
    hist.GetArray()[i] =this_norm > 2 ? val/this_norm : 0;
  }
}

inline double UCorrelator::Correlator::getDeltaT(int ant1, int ant2, double phi, double theta)
{
  double part1=ap->zAntByPolAnt[pol][ant1]*tan(theta * DEG2RAD) - ap->rAntByPolAnt[pol][ant1] * cos(phi *DEG2RAD);
  double part2=ap->zAntByPolAnt[pol][ant2]*tan(theta * DEG2RAD) - ap->rAntByPolAnt[pol][ant2] * cos(phi *DEG2RAD); 
  
  double geomDelay=1e9*((cos(theta *DEG2RAD ) * (part1 - part2))/C_LIGHT);    //returns time in ns
  
  return geomDelay;
}

inline double UCorrelator::Correlator::getDeltaTFast(int ant1, int ant2, int phibin, int thetabin)
{
  double part1=ap->zAntByPolAnt[pol][ant1]*tan_theta_table[thetabin] - ap->rAntByPolAnt[pol][ant1] * cos_phi_table[phibin];
  double part2=ap->zAntByPolAnt[pol][ant2]*tan_theta_table[thetabin] - ap->rAntByPolAnt[pol][ant2] * cos_phi_table[phibin];
  
  double geomDelay=1e9*((cos_theta_table[thetabin] * (part1 - part2))/C_LIGHT);    //returns time in ns

  return geomDelay;
}


inline Double_t UCorrelator::Correlator::getGroupDelay(Double_t phiToAntBoresight, Double_t thetaWave)//in radians, with positive being down
{
  double thetaDeg=thetaWave-10;
  double phiDeg=phiToAntBoresight;
  Double_t totalAngleDeg=sqrt(thetaDeg*thetaDeg+phiDeg*phiDeg);//TODO
  if (totalAngleDeg>50) totalAngleDeg=50;
  //Double_t delayTime=(totalAngleDeg*totalAngleDeg*totalAngleDeg*totalAngleDeg)*1.303e-8;
  //delayTime+=(totalAngleDeg*totalAngleDeg*totalAngleDeg)*6.544e-8;
  //delayTime-=(totalAngleDeg*totalAngleDeg)*4.770e-6;
  //delayTime+=totalAngleDeg*8.097e-4;

  Double_t delayTime=(totalAngleDeg*totalAngleDeg*totalAngleDeg*totalAngleDeg)*1.45676e-8;
  delayTime-=(totalAngleDeg*totalAngleDeg)*5.01452e-6;

  return delayTime;
}

UCorrelator::Correlator::~Correlator()
{
  delete phi_table; 
  delete theta_table; 
  delete cos_phi_table; 
  delete tan_theta_table; 
  delete cos_theta_table; 

  reset(); 

}
