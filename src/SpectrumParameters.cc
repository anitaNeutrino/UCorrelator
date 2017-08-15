#include "SpectrumParameters.h" 
#include "TGraph.h" 
#include "TLinearFitter.h"
#include "TMutex.h"
#include "AnalysisConfig.h"
#include <cfloat>
#include "TMath.h"
#ifdef UCORRELATOR_OPENMP
#include <omp.h>
#endif



/** This algorithm is kinda crappy right now. 
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 */ 

static TLinearFitter ** fitters; 

static void setupFitters() __attribute__((constructor));

void setupFitters() 
{
  int nthreads = 1; 

#ifdef UCORRELATOR_OPENMP
  nthreads = omp_get_max_threads(); 
#endif

//   printf("set up %d fitters\n", nthreads); 
  fitters = new TLinearFitter * [nthreads]; 
  for (int i = 0; i < nthreads; i++)
  {
    fitters[i] = new TLinearFitter(1,"1++x","x"); 
    double dummy_x[5] = {0}; 
    double dummy_y[5] = {0}; 
    /* !@#!DSAFSD #!@$ !@#!  %!@$ */ 
    fitters[i]->AssignData(5,1,dummy_x,dummy_y); 

  }
}

void UCorrelator::spectrum::fillSpectrumParameters(const TGraph * spectrum, const TGraph * average, 
                                                   AnitaEventSummary::WaveformInfo * winfo,
                                                   const AnalysisConfig * config) 
{
  TLinearFitter* fitter;

#ifdef UCORRELATOR_OPENMP
  fitter = fitters[omp_get_thread_num()]; 
#else
  fitter = fitters[0]; 
#endif


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


  //subtract off average spectrum 
  for (size_t i = 0; i < x.size(); i++) 
  {
    y[i] -= average->Eval(x[i]); 
  }

  double m = 0; 
  double b = 0; 

  std::vector<bool> used(x.size(), false); 

  for (int i = 0; i < AnitaEventSummary::peaksPerSpectrum; i++) 
  {
     double max_val = -DBL_MAX; 
     int max_j = -1; 

     for (unsigned j = 0; j < x.size(); j++) 
     {
       if (used[j]) continue; 
       double val =  y[j];// - (m * x[j] + b); 
//       printf("%d |   %f, %f\n", i, x[j], val); 
       if (val > max_val) 
       {
         max_j = j; 
         max_val = val; 
       }
     }

//     printf("\t max_i: %f\n", xx[max_i]); 

     int index_bounds[2] = {max_j,max_j}; 
     double power = TMath::Power(10,max_val/10); 

//      for (int sign = -1; sign <=1; sign+=2)
//      {
//         int how_far = 1; 
//         while(true) 
//         {
//           int jj = max_j + how_far * sign; 
//           double val = 0; 
// //          printf("%f\n",max_val - val); 
//           if (jj <= 0 || jj >=N-1 || max_val - val > config->bw_ndb || used[jj])
//           {
//             index_bounds[(sign+1)/2] = jj;
//             break; 
//           }
//           power += TMath::Power(10,val/10); 
//           how_far++; 
//         }
//      }
     winfo->peakPower[i] = 10 * log10(power); 

     int start = index_bounds[0];
     int end = index_bounds[1]; 
     if (start < 0) start = 0; 
     if (end >= (int) x.size()) end = x.size()-1; 

     winfo->bandwidth[i] = x[end] - x[start]; 
     winfo->peakFrequency[i] = (x[start] + x[end])/2;


//     printf("%d %d\n",start,end); 
//     printf("%f %f\n", power, x[max_j]); 
//     printf("%f %f\n", x[end]-x[start],winfo->peakFrequency[i] ); 

     
     //now remove adjacents up to ndb below 

     for (int u= start; u < end; u++) used[u] = true; 
     
   }


   // fit for slope 
   fitter->ClearPoints(); 

   //TODO do I need to mutex this?!? 
   fitter->AssignData(x.size(),1,&x[0] ,&y[0]); 

   fitter->Eval(); 
   m = fitter->GetParameter(1); 
   b = fitter->GetParameter(0); 

   winfo->spectrumSlope = m; 
   winfo->spectrumIntercept = b; 

}
