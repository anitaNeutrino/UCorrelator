#ifndef UCKDE_H_
#define UCKDE_H_
#include "TF2.h" 
#include "TH2.h" 
#include <vector> 

/** Kernel Density Estimator
 *  2d
 *
 **/ 

namespace UCorrelator
{
  class KDE2D 
  {

    public: 

      struct KDE2DOptions
      {
        //todo: implement ballooon 
        int adaptiveness; //0 for global, otherwise it uses the nth distance to scale sigma (sample-based) 
        double sigma_x;  //x width parameter, or 0 for auto 
        double sigma_y; //y width parameter, or 0 to auto
        double corr; //correlation (no auto setting yet, default 0); 


        KDE2DOptions(int adaptiveness = 0, double sigma_x = 0, double sigma_y = 0, double corr = 0) 
          : adaptiveness(adaptiveness), sigma_x(sigma_x), sigma_y(sigma_y), corr(corr) 
        { }
      };
      

      KDE2D(int n, const double * x, const double * y, const double * weights = 0,  KDE2DOptions = KDE2DOptions()); 

      double operator() (double x, double y) const; 
      double operator() (const double* X, const double * Par=0) const 
      {
        (void) Par; return this->operator()(X[0], X[1]); 
      }

      void getNSigmaBounds(double nsigma, double & xmin, double & xmax, double & ymin, double & ymax) const; 

      TF2 * makeTF2(const char * name = "kde2d", double nsigma_bound = 3) const; 
      /** make a histogram. by default, bin width will be auto_factor of smallest sigma x or y */ 
      TH2 * makeTH2(double binwidth_x = 0, double binwidth_y = 0,  const char * name = "hkde2d", const char * title = "Kernel Density Estimator", double n_sigma_bound = 3, double auto_factor = 0.25) const; 
      virtual ~KDE2D() { ;} 

      unsigned getN() const { return x.size(); } 
      double getX(int i) const { return x[i]; }
      double getY(int i) const { return y[i]; }
      double getWeight(int i) const { return y[i]; }

    private: 
      std::vector<double> x; 
      std::vector<double> y; 
      std::vector<double> w; 
      std::vector<double> sigma_x; 
      std::vector<double> sigma_y; 
      std::vector<double> rho; 
      double W; 

      ClassDef(KDE2D,1); 
  };



}

#endif

