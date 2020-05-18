#include "CachedFC.h" 
#include <iostream> 
#include "Math/ProbFuncMathCore.h" 

#if __cplusplus <= 199711L


UCorrelator::CachedFC::CachedFC(double cl, int max_sig, double max_bg, double d_bg)
  : fc(cl), interp(0); 
{
  std::cerr << "CachedFC requires a C++11 compiler" << std::endl; 
}

double UCorrelator::CachedFC::upperLimit(int nsig, double bg) 
{
  return -1; 

}

UCorrelator::CachedFC::~CachedFC() 
{
}

#else

#include "lazylookup.hh" 

UCorrelator::CachedFC::CachedFC(double cl, int max_sig, double max_bg, double d_bg)
  : fc(cl)
{

  //figure out what MuMax should be
  int mu_max = ceil(max_bg); 
  while (ROOT::Math::poisson_cdf_c(mu_max, max_bg) > 1e-9)
  {
    mu_max++; 
  }

  fc.SetMuMax(mu_max); 


  interp = (void*)  new lazylookup::grid_interpolator<2>( [&](const std::array<double,2>& X) { return fc.CalculateUpperLimit(X[0], X[1]); }, 
                                                 {(unsigned) max_sig, (unsigned) (max_bg/d_bg)}, {0.,0.},{(double) max_sig, max_bg}); 

}

double UCorrelator::CachedFC::upperLimit(int nsig, double bg) 
{
  double dsig= nsig; 
  return ( (lazylookup::grid_interpolator<2>*) interp)->val({dsig,bg}); 

}

UCorrelator::CachedFC::~CachedFC() 
{
  delete (lazylookup::grid_interpolator<2>*) interp; 
}
#endif
