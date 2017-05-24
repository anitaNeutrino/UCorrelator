#include "PointingResolutionModel.h" 
#include "TMath.h" 



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


double UCorrelator::PointingResolution::computeProbability(double _phi, double _theta)
{
  double ans; 
  computeProbability(1,&_phi,&_theta,&ans) ; 
  return ans; 
}



double * UCorrelator::PointingResolution::computeProbability(int N, 
                                                             const double * __restrict__ vphi,
                                                             const double * __restrict__ vtheta, 
                                                             double * __restrict__ p)
{
  if (!p) p = new double[N]; 

  for (int i = 0; i < N; i++) 
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
