#ifndef UCORRELATOR_RESOLUTION_MODEL_H
#define UCORRELATOR_RESOLUTION_MODEL_H

#include "AnitaEventSummary.h"

namespace UCorrelator
{

  class PointingResolution 
  {
    public: 

      PointingResolution(double phi = 0, double theta = 0, double dphi = 0, double dtheta = 0, double rho = 0); 
      virtual ~PointingResolution() {;} 

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
      ClassDef(PointingResolution,1); 

  }; 

  class PointingResolutionModel 
  {
    public: 
      virtual PointingResolution *  computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution * p = 0) const = 0; 
      virtual ~PointingResolutionModel() { ; } 
      ClassDef(PointingResolutionModel,1); 
  }; 



  class ConstantPointingResolutionModel : public PointingResolutionModel
  {
    public: 
      ConstantPointingResolutionModel(double dphi=0.4, double dtheta=0.3, double rho = 0) 
        : dphi(dphi), dtheta(dtheta), rho(rho) { }


      virtual PointingResolution * computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution *p) const
      {
        new (p) PointingResolution(sum->peak[pol][peak].phi, sum->peak[pol][peak].theta, dphi,dtheta,rho); 
        return p; 
      }

      private: 
        double dphi, dtheta, rho; 
        ClassDef(ConstantPointingResolutionModel,1); 

  }; 

}




#endif 
