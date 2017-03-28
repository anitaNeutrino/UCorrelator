#include "SystemResponse.h" 
#include "FFTtools.h" 
#include "TF1.h" 

#include "AnalysisWaveform.h" 


static UCorrelator::BandLimitedDeconvolution bld(0.18,1.1); 
UCorrelator::DeconvolutionMethod & UCorrelator::kDefaultDeconvolution = bld; 

void UCorrelator::WienerDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  FFTWComplex zero(0,0); 
  for (unsigned i = 0; i < N; i++) 
  {
    double f = i * df; 

    double SNR = snr(f); 
//    printf("SNR(%f) = %f\n", f, SNR); 

    if (SNR <=0)  Y[i] = zero; 
    else
    {
      double H2 = response[i].getAbsSq(); 
      Y[i] = Y[i] / response[i] * ( H2 / (H2 + 1./SNR)); 
    }
  }
}


UCorrelator::WienerDeconvolution::WienerDeconvolution(const TGraph *g, const double * sc) 
{
  snr_graph = g; 
  min = g->GetX()[0]; 
  max = g->GetX()[g->GetN()-1]; 
  snr_function = 0; 
  scale = sc; 
}

UCorrelator::WienerDeconvolution::WienerDeconvolution(const TF1 *f) 
{
  snr_function = f; 
  f->GetRange(min,max); 
  snr_graph = 0;
  scale = 0; 
}



double UCorrelator::WienerDeconvolution::snr(double f) const
{

  if (f < min || f > max ) return 0; 

  if (snr_graph)
  {
    double sc = scale ? *scale : 1; 
    return sc*snr_graph->Eval(f); 
  }
  else if (snr_function) 
  {
//    printf("%f\n",f); 
    return snr_function->Eval(f); 
  }
 
  fprintf(stderr,"Something's wrong...\n"); 

  return -1; 
}

void UCorrelator::AllPassDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  (void) df; 
  for (unsigned i = 0; i < N; i++) 
  {
    FFTWComplex r = response[i]; 
    r/= r.getAbs(); 
    Y[i]/=r; 
  }
}



void UCorrelator::NaiveDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  (void) df; 
  for (unsigned i = 0; i < N; i++) 
  {
    Y[i]/=response[i]; 
  }
}

void UCorrelator::BandLimitedDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  size_t min_i =(size_t) TMath::Max(0,int(min_freq / df)); 
  size_t max_i =(size_t) TMath::Min((int) N-1, int(max_freq / df)); 
  FFTWComplex zero(0,0); 

  for (size_t i = 0; i < min_i; i++) 
  {
    Y[i] = zero; 
  }

  for (size_t i = min_i; i < max_i; i++) 
  {
    Y[i]/=response[i]; 
  }
  for (size_t i = max_i; i < N; i++) 
  {
    Y[i] = zero; 
  }

}



UCorrelator::Response::Response(int Nfreq, double df, int nangles, const double * the_angles, const FFTWComplex ** the_responses)  
  : Nfreq(Nfreq), df(df)
{
  for (int i = 0; i < nangles; i++) 
  {
    addResponseAtAngle(the_angles[i], the_responses[i]); 
  }

}

void UCorrelator::Response::addResponseAtAngle(double angle, const FFTWComplex * response) 
{
    FFTWComplex * copy = new FFTWComplex[Nfreq]; 
    memcpy(copy, response, Nfreq * sizeof(FFTWComplex)); 
    responses[angle] = copy; 
//    if (angle == 0) response_at_zero = copy; 
    dirty = true; 
}

UCorrelator::Response::Response(int Nfreq, double df)
 : Nfreq(Nfreq), df(df) 
{

}

UCorrelator::Response::Response(const TGraph * timedomain, int npad)
{
  AnalysisWaveform aw(timedomain->GetN(), timedomain->GetX(), timedomain->GetY(), timedomain->GetX()[1] - timedomain->GetX()[0]); 
  aw.padEven(npad); 
  (void)  aw.freq(); 
  df = aw.deltaF(); 
  Nfreq = aw.Nfreq(); 
  addResponseAtAngle(0,aw.freq()); 
}



UCorrelator::Response::Response(int Nfreq, double df, const FFTWComplex * response)  
  : Nfreq(Nfreq),df(df)
{
  addResponseAtAngle(0, response); 
}


FFTWComplex * UCorrelator::AbstractResponse::getResponseArray(int N, const double * f, double angle) const
{
  FFTWComplex * answer = new FFTWComplex[N]; 
  for (int i = 0; i < N; i++) answer[i] = getResponse(f[i], angle); 
  return answer; 

}

FFTWComplex * UCorrelator::AbstractResponse::getResponseArray(int N, double  df, double angle) const
{
  FFTWComplex * answer = new FFTWComplex[N]; 
  for (int i = 0; i < N; i++) answer[i] = getResponse(i*df, angle); 
  return answer; 
}



void UCorrelator::Response::recompute()  const
{

  int nangles = responses.size(); 
  real.SetBins( Nfreq, 0, df * Nfreq,nangles, -90, 90); 
  imag.SetBins( Nfreq, 0, df * Nfreq, nangles, -90, 90); 
  

 
  if (nangles > 1) 
  {
    double bin_boundaries[nangles+1]; 
    double center_angles[nangles]; 

   //fill in centers of angles
    int i = 0; 
    for (std::map<double, FFTWComplex *>::const_iterator it = responses.begin(); it!=responses.end(); it++)
    {
      center_angles[i++] = it->first; 
    }

    
    bin_boundaries[0] = center_angles[0] - (center_angles[1]-center_angles[0])/2;    
    bin_boundaries[nangles] = center_angles[nangles-1] + (center_angles[nangles-1]-center_angles[nangles-2])/2;    

    for (int i = 1; i < nangles; i++) 
    {
      bin_boundaries[i] = (center_angles[i-1] + center_angles[i])/2; 
    }

    imag.GetYaxis()->Set(nangles, bin_boundaries); 
    real.GetYaxis()->Set(nangles, bin_boundaries); 
  }

  int j = 1; 

  for (std::map<double, FFTWComplex *>::const_iterator it = responses.begin(); it!=responses.end(); it++)
  {
    for (int i = 0; i <  Nfreq; i++) 
    {
      imag.SetBinContent(i+1,j, it->second[i].im); 
      real.SetBinContent(i+1,j, it->second[i].re); 
    }
    j++; 
  }

//  imag.Print(); 

  dirty = false; 
}

FFTWComplex UCorrelator::CompositeResponse::getResponse(double f, double angle )  const
{

  FFTWComplex answer(1,0); 
  for (size_t i = 0; i < responses.size(); i++) answer*= responses[i]->getResponse(f,angle); 
  return answer; 
}






FFTWComplex UCorrelator::Response::getResponse(double f, double angle ) const
{
  lock.Lock();
  if (dirty)
  {
    recompute(); 
  }
  lock.UnLock(); 

//  printf("%f %f %f %f %f %f\n", f, angle, real.GetXaxis()->GetXmin(), real.GetXaxis()->GetXmax(), real.GetYaxis()->GetXmin(), real.GetYaxis()->GetXmax()); 
  double re = real.Interpolate(f,angle); 
  double im = imag.Interpolate(f,angle); 
//  printf("%f %f %f\n", f, re, im); 
  
  return FFTWComplex(re,im); 
}

double UCorrelator::AbstractResponse::getMagnitude(double f, double angle) const 
{
  return getResponse(f,angle).getAbs(); 
}

double UCorrelator::AbstractResponse::getPhase(double f, double angle) const 
{
  return getResponse(f,angle).getPhase(); 
}






AnalysisWaveform* UCorrelator::AbstractResponse::impulseResponse(double dt, int N )  const
{
  AnalysisWaveform * out = new AnalysisWaveform(N, dt); 
  out->updateEven()->GetY()[1] = 1; 
  convolveInPlace(out,0); 
  return out; 
}

void UCorrelator::AbstractResponse::convolveInPlace(AnalysisWaveform * wf, double angle)  const
{
  int old_size = wf->Neven(); 
//  wf->padEven(2); 
  int nf = wf->Nfreq(); 
  double df = wf->deltaF(); 
  FFTWComplex * fft = wf->updateFreq(); 
  for (int i = 0; i < nf; i++) 
  {
    fft[i] *= getResponse(i*df, angle); 
  }
  wf->updateEven()->Set(old_size); 
}

AnalysisWaveform * UCorrelator::AbstractResponse::convolve(const AnalysisWaveform * in,  double off_axis_angle ) const
{

  //copy 
  AnalysisWaveform * out = new AnalysisWaveform(*in); 
  convolveInPlace(out, off_axis_angle); 
  return out; 
}

AnalysisWaveform * UCorrelator::AbstractResponse::deconvolve(const AnalysisWaveform * in,  const DeconvolutionMethod * method, double off_axis_angle) const
{

  //copy 
  AnalysisWaveform * out = new AnalysisWaveform(*in); 
  deconvolveInPlace(out,method,off_axis_angle); 
  return out; 
}


/** caching vars */ 
static __thread const UCorrelator::AbstractResponse *cache_response = 0; 
static __thread int cache_Nf= 0; 
static __thread double cache_df= 0; 
static __thread double cache_angle= 0; 
static __thread FFTWComplex *cache_V = 0; 


void UCorrelator::AbstractResponse::deconvolveInPlace(AnalysisWaveform * wf,  const DeconvolutionMethod * method, double off_axis_angle) const
{
//  printf("method: %p\n", method); 
  int old_size = wf->Neven(); 
  wf->padEven(1); 
  int nf = wf->Nfreq();
  double df = wf->deltaF(); 
  if (!cache_response || cache_Nf != nf || cache_df != df || !cache_V || cache_angle != off_axis_angle) 
  {
    cache_response = this; 
    cache_Nf = nf; 
    cache_df = df; 
    cache_angle = off_axis_angle; 
    if (cache_V) delete [] cache_V; 
    cache_V = getResponseArray(nf,df,off_axis_angle);
  }
    
  method->deconvolve(nf,df, wf->updateFreq(), cache_V); 

  wf->updateEven()->Set(old_size); 

}
