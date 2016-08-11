#include "SystemResponse.h" 
#include "FFTtools.h" 

#include "AnalysisWaveform.h" 




void UCorrelator::NaiveDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  for (unsigned i = 0; i < N; i++) 
  {
    Y[i]/=response[i]; 
  }
}

void UCorrelator::BandLimitedDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  int min_i = TMath::Max(0,int(min_freq / df)); 
  int max_i = TMath::Min(N-1, size_t(max_freq / df)); 
  for (int i = min_i; i < max_i; i++) 
  {
    Y[i]/=response[i]; 
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


/*
static FFTWComplex interpolateBetween(const FFTWComplex *a, const FFTWComplex * b, double frac) 
{

  double mag = a->getAbs() *(1- frac) + frac * b->getAbs();

  double phaseA = a->getPhase(); 
  double phaseB = b->getPhase(); 

  double dphase = phaseB-phaseA; 
  if (dphase > TMath::Pi()) dphase -= 2*TMath::Pi(); 
  if (dphase < -TMath::Pi()) dphase += 2*TMath::Pi(); 

  double phase = phaseA + dphase * frac; 

  return FFTWComplex(mag*cos(phase), mag*sin(phase)); 
}
 unneeded right now */


void UCorrelator::Response::recompute()  const
{

  int nangles = responses.size(); 
  phases.SetBins( Nfreq, 0, df * Nfreq,nangles, -90, 90); 
  mags.SetBins( Nfreq, 0, df * Nfreq, nangles, -90, 90); 
  
 
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

    phases.GetYaxis()->Set(nangles, bin_boundaries); 
    mags.GetYaxis()->Set(nangles, bin_boundaries); 
  }

  int j = 1; 

  for (std::map<double, FFTWComplex *>::const_iterator it = responses.begin(); it!=responses.end(); it++)
  {
    for (int i = 0; i <= Nfreq; i++) 
    {
      phases.SetBinContent(i+1,j, it->second[i].getPhase() ); 
      mags.SetBinContent(i+1,j, it->second[i].getAbs() ); 
    }
    j++; 
  }
}

FFTWComplex UCorrelator::CompositeResponse::getResponse(double f, double angle )  const
{

  FFTWComplex answer(1,0); 
  for (size_t i = 0; i < responses.size(); i++) answer*= responses[i]->getResponse(f,angle); 
  return answer; 
}




double UCorrelator::Response::getPhase(double f, double angle ) const
{
  if (dirty)
  {
    recompute(); 
  }

  return phases.Interpolate(f,angle); 
}


double UCorrelator::Response::getMagnitude(double f, double angle ) const
{
  if (dirty)
  {
    recompute(); 
  }

  return mags.Interpolate(f,angle);  
}



FFTWComplex UCorrelator::Response::getResponse(double f, double angle ) const
{

  double mag = getMagnitude(f,angle); 
  double phase = getPhase(f,angle); 
  return FFTWComplex(mag * cos(phase) , mag * sin(phase)); 
}

double UCorrelator::AbstractResponse::getMagnitude(double f, double angle) const 
{
  return getResponse(f,angle).getAbs(); 
}

double UCorrelator::AbstractResponse::getPhase(double f, double angle) const 
{
  return getResponse(f,angle).getPhase(); 
}







void UCorrelator::AbstractResponse::convolveInPlace(AnalysisWaveform * wf, double angle)  const
{
  int old_size = wf->Neven(); 
  wf->padEven(2); 
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

void UCorrelator::AbstractResponse::deconvolveInPlace(AnalysisWaveform * wf,  const DeconvolutionMethod * method, double off_axis_angle) const
{
  wf->padEven(2); 
  int nf = wf->Nfreq();
  double df = wf->deltaF(); 
  std::vector<FFTWComplex> R(nf); 
  for (int i = 0; i < nf; i++) 
  {
    R[i] =getResponse(i * df, off_axis_angle);
  }
    
  method->deconvolve(nf,df, wf->updateFreq(), &R[0]); 

}
