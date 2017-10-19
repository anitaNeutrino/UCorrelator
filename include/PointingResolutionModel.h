#ifndef UCORRELATOR_RESOLUTION_MODEL_H
#define UCORRELATOR_RESOLUTION_MODEL_H

#include "AnitaEventSummary.h"
class TRandom; 
class TF1; 

namespace UCorrelator
{

  class PointingResolution 
  {
    public: 

      PointingResolution(double phi = 0, double theta = 0, double dphi = 0, double dtheta = 0, double rho = 0); 
      virtual ~PointingResolution() {;} 

      double * computeProbabilityDensity(int N, const double *phi, 
                                         const double *theta,
                                         double * p = 0); 
      void random(double *phi, double * theta, TRandom * rng= 0);  

      double computeProbabilityDensity(double phi, double theta); 

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



  
  class PointingResolutionParSNRModel : public PointingResolutionModel
  {
    public: 
     PointingResolutionParSNRModel(TF1 * f_dtheta, TF1 * f_dphi, bool use_deconvolved = false)
     : f_th(f_dtheta), f_ph(f_dphi), deconv(use_deconvolved)  {; } 
       ; 

     virtual PointingResolution * computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution *p) const; 

    private: 
     TF1 * f_th; 
     TF1 * f_ph; 
     bool deconv; 
      ClassDef(PointingResolutionParSNRModel,1); 

  }; 


  class ConstantPointingResolutionModel : public PointingResolutionModel
  {
    public: 
      ConstantPointingResolutionModel(double dphi=0.4, double dtheta=0.3) 
        : dphi(dphi), dtheta(dtheta) { rho = 0; }


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
