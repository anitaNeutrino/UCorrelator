#ifndef UCORRELATOR_RESOLUTION_MODEL_H
#define UCORRELATOR_RESOLUTION_MODEL_H

#include "AnitaEventSummary.h"
#include "TF1.h" 
class TRandom; 
class TProfile; 

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
      double getdPhi() const { return dphi ; }
      double getdTheta() const { return dtheta ; }
      double getCorr() const { return rho; } 


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



  class HeadingErrorEstimator 
  {
    public: 
      HeadingErrorEstimator(int nseconds = 60)  : current_run(-1), nsecs(nseconds), prof(0) { ; } 
      int estimateHeadingError(double t,  double * stdev, double * offset = 0); 
      virtual ~HeadingErrorEstimator(); 

    private: 
      HeadingErrorEstimator(const HeadingErrorEstimator & );
      HeadingErrorEstimator operator= (const HeadingErrorEstimator &); 
      int current_run; 
      int nsecs; 
      TProfile * prof; 
  }; 

  class PointingResolutionModelPlusHeadingError : public PointingResolutionModel
  {
    public: 
     PointingResolutionModelPlusHeadingError(int nsecs, const PointingResolutionModel * other); 
     virtual PointingResolution * computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution *p) const; 
    private: 
     mutable HeadingErrorEstimator h; 
     const PointingResolutionModel *p; 
     ClassDef(PointingResolutionModelPlusHeadingError,1); 

  }; 
  
  class PointingResolutionParSNRModel : public PointingResolutionModel
  {
    public: 
     PointingResolutionParSNRModel() :  deconv(false){ ; } 
     PointingResolutionParSNRModel(const TF1 & f_dtheta, const TF1 & f_dphi, bool use_deconvolved = false, double scale_by_cos_theta = true)
     : f_th(f_dtheta), f_ph(f_dphi), deconv(use_deconvolved), cos_theta_scale(scale_by_cos_theta)  {; } 
       ; 

     virtual PointingResolution * computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution *p) const; 

    private: 
     TF1 f_th; 
     TF1 f_ph; 
     bool deconv; 
     bool cos_theta_scale; 
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
