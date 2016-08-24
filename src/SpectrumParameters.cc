#include "SpectrumParameters.h" 
#include "TGraph.h" 
#include "TLinearFitter.h"
#include "AnalysisConfig.h"
#include <cfloat>



/** This is a shitty first pass algorithm 
 *
 * Do line fit
 *  Find max subtracting
 *  Take out that peak and surrounding points 
 *  iterate 
 *
 *
 */ 

static __thread TLinearFitter * fitter = 0; 

void UCorrelator::spectrum::fillSpectrumParameters(const TGraph * spectrum, 
                                                   AnitaEventSummary::WaveformInfo * winfo,
                                                   const AnalysisConfig * config) 
{

  if (!fitter) fitter = new TLinearFitter(1,"1++x",""); 


  //TODO don't hardcode these 
  double lower_limit = config->spectral_fit_start; 
  double upper_limit = config->spectral_fit_stop; 

  double df = spectrum->GetX()[1]-spectrum->GetX()[0]; 
  double f0 = spectrum->GetX()[0]; 

  int low = (lower_limit - f0) / df; 
  int high = (upper_limit - f0) / df + 0.5; 
  if (high > spectrum->GetN()) high = spectrum->GetN(); 

  int N = high - low +1; 


  std::vector<double> x(N); 
  std::vector<double> y(N); 

  const double * xx = spectrum->GetX() + low; 
  const double * yy = spectrum->GetY() + low; 

  memcpy(&x[0], xx, N * sizeof(double)); 
  memcpy(&y[0], yy, N * sizeof(double)); 

  double m = 0; 
  double b = 0; 

  std::vector<int> maxes; 

  for (int i = 0; i < AnitaEventSummary::peaksPerSpectrum; i++) 
  {
    fitter->ClearPoints(); 
    fitter->AssignData(x.size(),1,&x[0] ,&y[0]); 
    fitter->Eval(); 
    m = fitter->GetParameter(1); 
    b = fitter->GetParameter(0); 

    double max_val = -DBL_MAX; 
    int max_i = -1; 

    for (unsigned j = 0; j < x.size(); j++) 
    {
      double val =  y[j] - (m * x[j] + b); 

      if (val > max_val) 
      {
        max_val = val; 
        max_i = j; 
      }
    }

    maxes.push_back(max_i); 

    //now remove adjacents 

    int start = max_i > 0 ? max_i-1 : 0; 
    int end = max_i+1; 

    x.erase(x.begin() + start, x.begin() + end); 
    y.erase(y.begin() + start, y.begin() + end); 
  }

  //one last fit 
  fitter->ClearPoints(); 
  fitter->AssignData(x.size(),1,&x[0] ,&y[0]); 
  fitter->Eval(); 
  m = fitter->GetParameter(1); 
  b = fitter->GetParameter(0); 

  //now let's evalaute things 


  for (int i = 0; i < AnitaEventSummary::peaksPerSpectrum; i++) 
  {
    int j = maxes[i]; 
    winfo->peakFrequency[i] = xx[j]; 
    double max_val = yy[j] - (m*xx[j] + b); 
    int how_far = 1; 
    int index_bounds[2] = {j,j}; 
    double power = max_val; 

    for (int sign = -1; sign <=1; sign+=2)
    {
       int jj = j + how_far * sign; 
       if (jj < 0 || jj >= N) continue; 
       double val = yy[j] - (m*xx[jj]+ b); 
       if (max_val - val > config->bw_ndb)
       {
         index_bounds[(sign+1)/2] = jj;
         continue; 
       }
       power += val; 
       how_far++; 
    }

    winfo->peakPower[i] = power; 
    winfo->bandwidth[i] = xx[index_bounds[1]] - xx[index_bounds[0]]; 
  }

  winfo->spectrumSlope = m; 
  winfo->spectrumIntercept = b; 

}
