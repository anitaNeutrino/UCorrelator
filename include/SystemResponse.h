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

      virtual AnalysisWaveform * convolve(const AnalysisWaveform * wf, double angle = 0) ; 
      virtual AnalysisWaveform * deconvolve(const AnalysisWaveform * wf, double angle = 0) ; 
      virtual void convolveInPlace(AnalysisWaveform * wf, double angle = 0) ; 
      virtual void deconvolveinPlace(AnalysisWaveform * wf, double angle = 0) ; 

  };

  class Response : public AbstractResponse
  {
    public: 
      Response(int NFreq, double df); 
      Reponse(int Nfreq, double df, int nangles, const double * angles, const FFTWComplex ** responses);  
      Reponse(int Nfreq, double df, const FFTWComplex * response);  

      void addResponseAtAngle(double angle, const FFTWcomplex * response); 

      virtual FFTWComplex getResponse(double f, double angle = 0) const; 
      virtual double getMagnitude(double f, double angle= 0) const; 
      virtual double getPhase(double f, double angle = 0) const; 
       

      virtual AnalysisWaveform * convolve(const AnalysisWaveform * wf, double angle = 0) ; 
      virtual AnalysisWaveform * deconvolve(const AnalysisWaveform * wf, double angle = 0) ; 
      virtual void convolveInPlace(AnalysisWaveform * wf, double angle = 0) ; 
      virtual void deconvolveinPlace(AnalysisWaveform * wf, double angle = 0) ; 
      
      virtual ~Response(); 

    protected: 
      int Nfreq; 
      double df; 
      int nangles; 
      std::map<double, FFTWComplex *> responses; 
      TH2D phases; 
      TH2D mags; 
      bool dirty; 
      void recompute(); 
  }; 


  class CompositeResponse : public AbstractResponse
  {
    public:  
      void addResponse(const AbstractResponse * response) { response.push_back(response); } 
      virtual FFTWComplex getResponse(double f, double angle = 0) const; 
   
    private: 
      std::vector<const AbstractResponse * > responses; 
  }; 


  AnalysisWaveform * convolve(const AnalysisWaveform * wf, const AbstractResponse * response, double off_axis_angle); 
  AnalysisWaveform * deconvolve(const AnalysisWaveform * wf, const AbstractResponse * response,  const DeconvolutionMethod * method = &kDefaultDeconvolution, double off_axis_angle = 0); 



}



#endif 

