#ifndef _UCORRELATOR_SYSTEM_RESPONSE_H 
#define _UCORRELATOR_SYSTEM_RESPONSE_H

  /* oh boy, down the OO rabbit hole we go! 
   *
   *  this file defines methods for deconvolution and system responses
   *
   **/

#include <vector> 
#include <complex>
#include "TH2.h" 
#include "FFTWComplex.h" 
#include <map>
class AnalysisWaveform; 

namespace UCorrelator
{
  

 
  class DeconvolutionMethod
  {
  
    public: 
      virtual void deconvolve(size_t N, double df, FFTWComplex * Y, 
                              const FFTWComplex * response) const = 0; 

      virtual ~DeconvolutionMethod()  { ; } 
  }; 


  class NaiveDeconvolution : public DeconvolutionMethod
  {

    public: 

      virtual void deconvolve(size_t N, double df, FFTWComplex * Y, 
                              const FFTWComplex * response) const ; 

      virtual ~NaiveDeconvolution()  { ; } 

  }; 

  class BandLimitedDeconvolution : public DeconvolutionMethod
  {
    public: 

      BandLimitedDeconvolution(double min_freq, double max_freq, int edge_order = 0) 
        : min_freq(min_freq),max_freq(max_freq), edge_order(edge_order) 
      {
      }

      virtual void deconvolve(size_t N, double df, FFTWComplex * Y, 
                              const FFTWComplex * response) const; 

      virtual ~BandLimitedDeconvolution() { ; } 


    private: 
      double min_freq, max_freq;
      int edge_order;
  }; 


  static BandLimitedDeconvolution kDefaultDeconvolution(0.2,1.2); 



  class AbstractResponse
  {

   public: 
      virtual FFTWComplex getResponse(double f, double angle = 0) const = 0; 
      virtual FFTWComplex * getResponseArray(int N, const double  * f, double angle = 0) const = 0; 
      virtual FFTWComplex * getResponseArray(int N, double df, double angle = 0) const = 0; 
      virtual double getMagnitude(double f, double angle= 0) const;  
      virtual double getPhase(double f, double angle = 0) const; 

      virtual AnalysisWaveform * convolve(const AnalysisWaveform * wf, double angle = 0) const; 
      virtual AnalysisWaveform * deconvolve(const AnalysisWaveform * wf, const DeconvolutionMethod * method = &kDefaultDeconvolution, double angle = 0) const; 
      virtual void convolveInPlace(AnalysisWaveform * wf, double angle = 0) const; 
      virtual void deconvolveInPlace(AnalysisWaveform * wf, const DeconvolutionMethod * method = &kDefaultDeconvolution, double angle = 0) const; 

  };

  class Response : public AbstractResponse
  {
    public: 
      Response(int NFreq, double df); 
      Response(int Nfreq, double df, int nangles, const double * angles, const FFTWComplex ** responses);  
      Response(int Nfreq, double df, const FFTWComplex * response);  

      void addResponseAtAngle(double angle, const FFTWComplex * response); 

      virtual FFTWComplex getResponse(double f, double angle = 0) const; 
      virtual double getMagnitude(double f, double angle= 0) const; 
      virtual double getPhase(double f, double angle = 0) const; 
       
      
      virtual ~Response() { ; } 

    protected: 
      int Nfreq; 
      double df; 
      int nangles; 
      std::map<double, FFTWComplex *> responses; 
      mutable TH2D phases; 
      mutable TH2D mags; 
      mutable bool dirty; 
      void recompute() const; 
  }; 


  class CompositeResponse : public AbstractResponse
  {
    public:  
      void addResponse(const AbstractResponse * response) { responses.push_back(response); } 
      virtual FFTWComplex getResponse(double f, double angle = 0) const; 
      virtual ~CompositeResponse() { ; } 
   
    private: 
      std::vector<const AbstractResponse * > responses; 
  }; 


}



#endif 

