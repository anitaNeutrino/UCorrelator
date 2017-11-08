#include "PointingResolutionModel.h" 
#include "TMath.h" 
#include "TRandom.h" 
#include "TF1.h" 

#include "AnitaDataset.h" 
#include "TFile.h" 
#include "TTree.h" 
#include "TProfile.h" 

ClassImp(UCorrelator::PointingResolution); 
ClassImp(UCorrelator::PointingResolutionModel); 
ClassImp(UCorrelator::PointingResolutionModelPlusHeadingError); 
ClassImp(UCorrelator::ConstantPointingResolutionModel); 
ClassImp(UCorrelator::PointingResolutionParSNRModel); 

UCorrelator::PointingResolution::PointingResolution(double phi, double theta,
                                                    double dphi, double dtheta, double rho) 
  : phi(phi), theta(theta), dphi(dphi), dtheta(dtheta), rho(rho) 
{
  inv_dphi2 = 1./(dphi*dphi); 
  inv_dtheta2 = 1./ (dtheta*dtheta); 
  inv_dphidtheta = 2*rho/(dphi*dtheta); 
  expterm = -1./ ( 2 * (1-rho*rho)); 
  norm = 1./ ( 2 * TMath::Pi() * dphi * dtheta * sqrt ( 1-rho*rho)); 

}



double UCorrelator::PointingResolution::computeProbabilityDensity(double _phi, double _theta)
{
  double ans; 
  computeProbabilityDensity(1,&_phi,&_theta,&ans) ; 
  return ans; 
}



void UCorrelator::PointingResolution::random(double * p_phi, double * p_theta, TRandom * rng) 
{

  if (!rng) rng = gRandom; 

  double X = rng->Gaus(); 
  double Y = rng->Gaus(); 
  if (p_phi) *p_phi = phi + dphi * X; 
  if (p_theta) *p_theta = theta+dtheta * (rho * X + (1-rho) * Y); 
}


double * UCorrelator::PointingResolution::computeProbabilityDensity(int N, 
                                                             const double * __restrict__ vphi,
                                                             const double * __restrict__ vtheta, 
                                                             double * __restrict__ p)
{
  if (!p) p = new double[N]; 

  //TODO vectorize the first part of this loop
  
  int start = 0; 

#pragma omp simd 
  for (int i = start; i < N; i++) 
  {
    double phidiff = vphi[i]-phi; 
    double thetadiff = vtheta[i]-theta; 
    double z=  ( phidiff * phidiff * inv_dphi2) - 
               (phidiff * thetadiff * inv_dphidtheta) + 
               (thetadiff*thetadiff * inv_dtheta2); 
    z *= expterm; 
    p[i] = norm * exp(z); 
  }

  return p; 
}


UCorrelator::PointingResolution * UCorrelator::PointingResolutionParSNRModel::computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution *p) const 
{
  const AnitaEventSummary::WaveformInfo * info = deconv ? &sum->deconvolved[pol][peak] : &sum->coherent[pol][peak]; 
  double snr = info->snr; 
  double snr_min, snr_max; 
  f_ph.GetRange(snr_min, snr_max); 
  if (snr < snr_min) snr = snr_min; 
  if (snr > snr_max) snr = snr_max; 
  f_th.GetRange(snr_min, snr_max); 
  if (snr < snr_min) snr = snr_min; 
  if (snr > snr_max) snr = snr_max; 

  double scale = 1; 
  if (cos_theta_scale) 
  {
    scale = 1./cos(TMath::DegToRad() *  sum->peak[pol][peak].theta); 
  }

  new (p) PointingResolution(sum->peak[pol][peak].phi, sum->peak[pol][peak].theta, f_ph.Eval(snr) * scale, f_th.Eval(snr) * scale *scale, 0); 
  return p; 

}

static __thread int prof_counter = 0; 

int UCorrelator::HeadingErrorEstimator::estimateHeadingError(double t, double * stdev, double * offset) 
{
  int run = AnitaDataset::getRunAtTime(t); 

  if (current_run != run) 
  {
    TString str; 
    str.Form("%s/run%d/gpsEvent%d.root", AnitaDataset::getDataDir(), run,run); 
    TFile f(str); 
    if (!f.IsOpen()) return -1; 
    TTree * tree = (TTree*) f.Get("adu5PatTree"); 
    double min = tree->GetMinimum("realTime"); 
    double max = tree->GetMaximum("realTime"); 
    if (!prof) 
    {
      prof = new  TProfile(TString::Format("heading_prof_%d", prof_counter++), "heading profile", (max-min)/nsecs+2, min-nsecs, max+nsecs); 
      prof->SetErrorOption("s"); 
    }
    else
    {
      prof->Clear(); 
      prof->SetBins((max-min)/nsecs+2, min-nsecs, max+nsecs); 
      prof->SetDirectory(&f); 
    }
    tree->Draw(TString::Format("headingA-headingB:realTime>>%s", prof->GetName()), "weightA && weightB","profs goff");  
    prof->SetDirectory(0); 
    current_run = run; 
  }

  int bin = prof->GetXaxis()->FindBin(t); 
  if (stdev) *stdev = prof->GetBinError(bin); 
  if (offset) *offset = prof->GetBinContent(bin); 
  return prof->GetBinEntries(bin); 
}

UCorrelator::HeadingErrorEstimator::~HeadingErrorEstimator() 
{
  delete prof; 
}

UCorrelator::PointingResolutionModelPlusHeadingError::PointingResolutionModelPlusHeadingError(int nsecs, const PointingResolutionModel * other)
  : h(nsecs), p(other)
{


}

UCorrelator::PointingResolution * UCorrelator::PointingResolutionModelPlusHeadingError::computePointingResolution(const AnitaEventSummary * sum, AnitaPol::AnitaPol_t pol, int peak, PointingResolution * point) const
{
  point = p->computePointingResolution(sum,pol,peak,point); 

  double stdev, mean;
  int n = h.estimateHeadingError(sum->realTime, &stdev,&mean); 

  double dphi = point->getdPhi();
  if (!n) //no nearby GPS, Indiscriminately add 1 degree to dphi 
  {
    dphi = sqrt(dphi*dphi + 1); 
  }
  else
  {
    dphi = sqrt(dphi*dphi + stdev*stdev + mean*mean); //worsen the resolution 
  }

  new (point)  PointingResolution(sum->peak[pol][peak].phi, sum->peak[pol][peak].theta, dphi, point->getdTheta(), point->getCorr()); 

  return point; 
}


