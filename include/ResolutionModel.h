#ifndef UCORRELATOR_RESOLUTION_MODEL_H
#define UCORRELATOR_RESOLUTION_MODEL_H

#include "AnitaEventSummary.h"

namespace UCorrelator
{

  class PointingResolution 
  {
    public: 

      PointingResolution(double phi, double theta, double dphi, double dtheta, double rho = 0); 

      double * computeProbability(int N, const double *phi, 
                                         const double *theta,
                                         double * p = 0); 

      double computeProbability(double phi, double theta); 

    private:
      double phi; 
      double theta; 
      double dphi; 
      double dtheta; 
      double inv_dphi2; 
      double inv_dtheta2; 
      double inv_dphidtheta; 
      double rho; 
      double expterm; 
      double norm; 

  }; 

  class PointingResolutionModel 
  {
    virtual void computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution * p) = 0; 
  }; 


}




#endif 
