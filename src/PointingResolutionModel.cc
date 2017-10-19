#include "PointingResolutionModel.h" 
#include "TMath.h" 
#include "TRandom.h" 
#include "TF1.h" 


ClassImp(UCorrelator::PointingResolution); 
ClassImp(UCorrelator::PointingResolutionModel); 
ClassImp(UCorrelator::ConstantPointingResolutionModel); 
ClassImp(UCorrelator::PointingResolutionParSNRModel); 

UCorrelator::PointingResolution::PointingResolution(double phi, double theta,
                                                    double dphi, double dtheta, double rho) 
  : phi(phi), theta(theta), dphi(dphi), dtheta(dtheta), rho(rho) 
{
  inv_dphi2 = 1./(dphi*dphi); 
  inv_dtheta2 = 1./ (dtheta*dtheta); 
  inv_dphidtheta = 2*rho/(dphi*dtheta); 
  expterm = -1./ ( 2 * (1-rho*rho)); 
  norm = 1./ ( 2 * TMath::Pi() * dphi * dtheta * sqrt ( 1-rho*rho)); 

}


double UCorrelator::PointingResolution::computeProbabilityDensity(double _phi, double _theta)
{
  double ans; 
  computeProbabilityDensity(1,&_phi,&_theta,&ans) ; 
  return ans; 
}



void UCorrelator::PointingResolution::random(double * p_phi, double * p_theta, TRandom * rng) 
{

  if (!rng) rng = gRandom; 

  double X = rng->Gaus(); 
  double Y = rng->Gaus(); 
  if (p_phi) *p_phi = phi + dphi * X; 
  if (p_theta) *p_theta = theta+dtheta * (rho * X + (1-rho) * Y); 
}


double * UCorrelator::PointingResolution::computeProbabilityDensity(int N, 
                                                             const double * __restrict__ vphi,
                                                             const double * __restrict__ vtheta, 
                                                             double * __restrict__ p)
{
  if (!p) p = new double[N]; 

  //TODO vectorize the first part of this loop
  
  int start = 0; 

#pragma omp simd 
  for (int i = start; i < N; i++) 
  {
    double phidiff = vphi[i]-phi; 
    double thetadiff = vtheta[i]-theta; 
    double z=  ( phidiff * phidiff * inv_dphi2) - 
               (phidiff * thetadiff * inv_dphidtheta) + 
               (thetadiff*thetadiff * inv_dtheta2); 
    z *= expterm; 
    p[i] = norm * exp(z); 
  }

  return p; 
}


UCorrelator::PointingResolution * UCorrelator::PointingResolutionParSNRModel::computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution *p) const 
{
  const AnitaEventSummary::WaveformInfo * info = deconv ? &sum->deconvolved[pol][peak] : &sum->coherent[pol][peak]; 
  double snr = info->snr; 
  double snr_min, snr_max; 
  f_ph->GetRange(snr_min, snr_max); 
  if (snr < snr_min) snr = snr_min; 
  if (snr < snr_max) snr = snr_max; 
  f_th->GetRange(snr_min, snr_max); 
  if (snr < snr_min) snr = snr_min; 
  if (snr < snr_max) snr = snr_max; 

  new (p) PointingResolution(sum->peak[pol][peak].phi, sum->peak[pol][peak].theta, f_ph->Eval(snr), f_th->Eval(snr), 0); 
  return p; 

}
