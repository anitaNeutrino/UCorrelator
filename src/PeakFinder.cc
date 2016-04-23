#include "PeakFinder.h"
#include "TH2.h" 
#include <Eigen/Dense> 
#include "TF2.h" 
#include "FFTtools.h"
#include "Math/Interpolator.h"
#include "TGraph.h"


static bool isLocalMaxima(const TH2D * hist, int bin) 
{


  const double* I = hist->GetArray(); 
  const int width = hist->GetNbinsX()+2; 

  // in the unlikely case that there are two bins with exactly the same value 
  // we'll pick the one to the top left because why not
  // I haven't checked if that really fixes it but oh well 
  if (I[bin-1] >= I[bin]) return false; 
  if (I[bin+1] > I[bin]) return false; 
  if (I[bin-width-1] >= I[bin]) return false;
  if (I[bin-width] >= I[bin]) return false;
  if (I[bin-width+1] >= I[bin]) return false;
  if (I[bin+width+1] > I[bin]) return false;
  if (I[bin+width] > I[bin]) return false;
  if (I[bin+width+1] > I[bin]) return false;


  return true; 
}

static int findMaximum(const TH2D* hist, std::vector<bool> * notallowed)
{
  const double * arr = hist->GetArray(); 
  int width = hist->GetNbinsX() + 2; 
  int height = hist->GetNbinsY() + 2; 
  int N = width * height;

  double max = -DBL_MAX; 
  int max_i = -1; 
  for (int i = 0; i < N; i++) 
  {
    if (notallowed->at(i))
    {
 //     printf("Skipping bin %d\n",i); 
      continue; 
    }
    //skip overflow
    
    if (i % width == 0) continue; 
    if (i % width == width-1) continue; 
    if (i / width == 0) continue; 
    if (i / width == height-1) break; 

    if (arr[i] > max && isLocalMaxima(hist,i))
    {
      max = arr[i]; 
      max_i = i; 
    }
  }
  return max_i; 
}

static void maskNearbyBins(const TH2D * hist, double distance, int bin, std::vector<bool> * used)
{

  const int width = hist->GetNbinsX()+2; 
  const int height = hist->GetNbinsY()+2; 
  int int_dist = int(distance+0.5); 

  int i0 = bin % width;
  int j0 = bin / width; 
  int dist2 = int(distance*distance+0.5);  

  for (int jj = -int_dist; jj <= int_dist; jj++)
  {
    int j = j0 +jj; 
    if (j < 1) continue; 
    if (j >= height) continue; 
    int j2 = jj*jj; 
    for (int ii = -int_dist; ii <= int_dist; ii++) 
    {
      int i = i0 +ii; 
      if (i < 1) continue; 
      if (i >= width) continue; 
      if (ii*ii + j2 > dist2) continue; 

      int bin2use =  i+ j * width; 
//      printf("Disallowing %d\n", bin2use); 
      (*used)[bin2use] = true; 
    }
  }
}


int UCorrelator::peakfinder::findIsolatedMaxima(const TH2D * hist, double distance, int Nmaxima, RoughMaximum * maxima,  bool use_bin_center)
{
  std::vector<bool> used( (hist->GetNbinsX() +2) * (hist->GetNbinsY()+2), false); 

  int nfound = 0; 

  while (nfound < Nmaxima)
  {
    int bin = findMaximum(hist, &used); 
    if (bin < 0) break; 

    maxima[nfound].bin = bin; 
    maxima[nfound].val = hist->GetArray()[bin]; 
    int xbin =  bin % (hist->GetNbinsX() + 2); 
    maxima[nfound].x = use_bin_center ? hist->GetXaxis()->GetBinCenter(xbin) : hist->GetXaxis()->GetBinLowEdge(xbin); 
    int ybin =  bin / (hist->GetNbinsX() + 2); 
    maxima[nfound].y = use_bin_center ? hist->GetYaxis()->GetBinCenter(ybin) : hist->GetYaxis()->GetBinLowEdge(ybin); 
    maskNearbyBins(hist, distance, bin, &used); 
    nfound++; 
  }

  for (int i = nfound; i < Nmaxima; i++)
  {
    maxima[i].bin  = -1; 
    maxima[i].val  = -9999; 
    maxima[i].x  = -9999; 
    maxima[i].y  = -9999; 
  }

  return nfound; 
}



void UCorrelator::peakfinder::FineMaximum::copyToPointingHypothesis(AnitaEventSummary::PointingHypothesis * p)
{
  p->phi = x; 
  p->theta = -y; 
  p->sigma_theta = sigma_y; 
  p->sigma_phi = sigma_x; 
  p->rho = covar / (sigma_x * sigma_y); 
  p->value = val; 
}

const int X2COEFF = 0; 
const int Y2COEFF = 1; 
const int XYCOEFF = 2; 
const int XCOEFF = 3; 
const int YCOEFF = 4; 
const int CCOEFF = 5; 


template <unsigned N> 
static const Eigen::JacobiSVD<Eigen::Matrix<double, N*N,6> > & getSVD()
{
  int halfway = N/2; 
  static Eigen::Matrix<double, N*N,6> M; 
  for (unsigned i = 0; i < N; i++)
  {
    for (unsigned j = 0; j < N; j++) 
    {
      int x =  i- halfway; 
      int y =  j - halfway; 
      M(i + j * N, X2COEFF) = x * x; 
      M(i + j * N, Y2COEFF) = y * y; 
      M(i + j * N, XYCOEFF) = x * y; 
      M(i + j * N, XCOEFF) = x; 
      M(i + j * N, YCOEFF) = y; 
      M(i + j * N, CCOEFF) = 1; 
    }
  }

  static Eigen::JacobiSVD<Eigen::Matrix<double, N*N,6> > svd(M,Eigen::ComputeFullU | Eigen::ComputeFullV); 
  return svd; 
}


template <unsigned N> 
static void doQuadraticPeakFinding(const TH2D * hist, UCorrelator::peakfinder::FineMaximum * peak)
{
  int xmax, ymax, ignored; 
  hist->GetMaximumBin(xmax,ymax,ignored); 

  double dx = hist->GetXaxis()->GetBinWidth(1); 
  double dy = hist->GetYaxis()->GetBinWidth(1); 
  int nbinsx = hist->GetNbinsX(); 
  int nbinsy = hist->GetNbinsY(); 
  int width = nbinsx + 2; 

  //build the linear system
  static const Eigen::JacobiSVD<Eigen::Matrix<double, N*N, 6> > & svd = getSVD<N>();  //I  hope this works! 
  Eigen::Matrix<double,N*N,1> Z; 

  double xcenter = hist->GetXaxis()->GetBinCenter(xmax); 
  double ycenter = hist->GetYaxis()->GetBinCenter(ymax); 
  

  int halfway = N/2; 
  for (unsigned i = 0; i < N; i++)
  {
    for (unsigned j = 0; j < N; j++) 
    {
      int ii = xmax + i - halfway; 
      int jj = ymax + j - halfway; 
      double z = ii  < 1 || ii > nbinsx || jj < 1 || jj > nbinsy ? 0 : hist->GetArray()[ii + jj * width]; 
      Z(i + j * N)  = z; 
    }
  }

  //solve the linear system
  Eigen::Matrix<double,6,1> B = svd.solve(Z); 
    
  //now do some algebra to figure out what x0,y0,sigma_x,sigma_y, covar are 
  //someone needs to check this 
  double a = B(X2COEFF); 
  double b = B(Y2COEFF); 
  double c = B(XYCOEFF); 
  double d = B(XCOEFF); 
  double e = B(YCOEFF); 
  double f = B(CCOEFF); 

//  printf("[%f,%f,%f,%f,%f,%f]\n",a,b,c,d,e,f); 

  //from taking gradient of ax^2 + by^2 + cxy + dx + ey + f and setting equal to 0 
  double x0 = (2 * b * d / c - e) / (c - 4*a*b / c); 
  double y0 = ( -d - 2* a * x0)/ c; 

  peak->x = xcenter + x0*dx; 
  peak->y = ycenter + y0*dy; 

  // now use -hessian matrix to compute variance, covariance. Not sure this is right. 
  // Could instead try to match coefficients to taylor expansion, but that's a lot of work!
  peak->val =  a*x0*x0 + b * y0*y0 + c * x0*y0 + d * x0 + e * y0 + f; 
  peak->sigma_x =  sqrt(-2*b*peak->val / (4*a*b - c*c))*dx;
  peak->sigma_y =  sqrt(-2*a*peak->val / (4*a*b - c*c))*dy;
  peak->covar = c*peak->val / (4*a*b - c*c)*dx*dy; 

}


void UCorrelator::peakfinder::doPeakFindingQuadratic9(const TH2D * hist, FineMaximum * peak) 
{
  doQuadraticPeakFinding<3>(hist,peak); 
}

void UCorrelator::peakfinder::doPeakFindingQuadratic16(const TH2D * hist, FineMaximum * peak) 
{
  doQuadraticPeakFinding<4>(hist,peak); 
}

void UCorrelator::peakfinder::doPeakFindingQuadratic25(const TH2D * hist, FineMaximum * peak) 
{
  doQuadraticPeakFinding<5>(hist,peak); 
}

void UCorrelator::peakfinder::doPeakFindingQuadratic36(const TH2D * hist, FineMaximum * peak) 
{
  doQuadraticPeakFinding<6>(hist,peak); 
}

void UCorrelator::peakfinder::doPeakFindingQuadratic49(const TH2D * hist, FineMaximum * peak) 
{
  doQuadraticPeakFinding<7>(hist,peak); 
}



void UCorrelator::peakfinder::doInterpolationPeakFindingAbby(const TH2D * hist, FineMaximum * peak)
{
  // This is mostly Abby's code, with some stuff in front and after to make it work here
  // As far as I can tell, if there are any correlations between theta and phi in the peak, this will not give you 
  // the right answer, but maybe I'm wrong. 
  
  int peakPhiBin, peakThetaBin, ignored; 
  int max_bin = hist->GetMaximumBin(peakPhiBin,peakThetaBin,ignored); 
  double peakVal = hist->GetArray()[max_bin]; 

  int NUM_BINS_FINE_THETA = hist->GetNbinsY(); 
  int NUM_BINS_FINE_PHI = hist->GetNbinsX(); 
  double FINE_BIN_SIZE=hist->GetXaxis()->GetBinWidth(1); 

  int npointsTheta=NUM_BINS_FINE_THETA - 19; // ok... 
  int npointsPhi=NUM_BINS_FINE_THETA - 19;// why not... 

  double thetaArray[npointsTheta];
  double phiArray[npointsPhi];
  double thetaValArray[npointsTheta], phiValArray[npointsPhi];

  for (int i=0;i<npointsTheta;i++){//construct an array of 10 values in phi and theta around the peak
    if (int(peakThetaBin-npointsTheta/2.)>=0 && int(peakThetaBin+npointsTheta/2.)<NUM_BINS_FINE_THETA){
      thetaValArray[i]=hist->GetBinContent(int(peakThetaBin-npointsTheta/2+i),peakPhiBin);
      thetaArray[i]=hist->GetYaxis()->GetBinCenter(int(peakThetaBin-npointsTheta/2+i));
    }
    
    else if(int(peakThetaBin-npointsTheta/2.)<0){
      thetaValArray[i]=hist->GetBinContent(i,peakPhiBin);
      thetaArray[i]=hist->GetYaxis()->GetBinCenter(i);
    }
    else{
      thetaValArray[i]=hist->GetBinContent(NUM_BINS_FINE_THETA-(npointsTheta-i),peakPhiBin);
      thetaArray[i]=hist->GetYaxis()->GetBinCenter(NUM_BINS_FINE_THETA-(npointsTheta-i)); 
    }
  }    
  for (int i=0;i<npointsPhi;i++){
    if (int(peakPhiBin-npointsPhi/2.)>=0 && int(peakPhiBin+npointsPhi/2.)<NUM_BINS_FINE_PHI){
      phiValArray[i]=hist->GetBinContent(peakThetaBin,int(peakPhiBin-npointsPhi/2+i));
      phiArray[i]=hist->GetXaxis()->GetBinCenter(int(peakPhiBin-npointsPhi/2+i));
    }
    
    else if (int(peakPhiBin-npointsPhi/2.)<0){ 
      phiValArray[i]=hist->GetBinContent(peakThetaBin,i);
      phiArray[i]=hist->GetXaxis()->GetBinCenter(i); 
    }
    else{
      phiValArray[i]=hist->GetBinContent(peakThetaBin,NUM_BINS_FINE_PHI-(npointsPhi-i));
      phiArray[i]=hist->GetXaxis()->GetBinCenter(NUM_BINS_FINE_PHI-(npointsPhi-i));
    }
    //cout<<"i: "<<i<<", phi: "<<phiArray[i]<<", phiVal: "<<phiValArray[i]<<endl;
    //cout<<"i: "<<i<<", theta: "<<thetaArray[i]<<", thetaVal: "<<thetaValArray[i]<<endl;
  }

  //do interpolation with akima to get peak and find FWHM or something like that.
  std::vector<double> thetaVector(npointsTheta);
  std::vector<double> thetaValVector(npointsTheta);
  std::vector<double> phiVector(npointsPhi);
  std::vector<double> phiValVector(npointsPhi);
  for (int i=0;i<npointsTheta;i++){
    thetaVector[i]=thetaArray[i];
    thetaValVector[i]=thetaValArray[i];
  }
 for (int i=0;i<npointsPhi;i++){
   phiVector[i]=phiArray[i];
   phiValVector[i]=phiValArray[i];
 }

  ROOT::Math::Interpolator interpolatorTheta(thetaVector.size(), ROOT::Math::Interpolation::kAKIMA );
  interpolatorTheta.SetData(thetaVector,thetaValVector);
  ROOT::Math::Interpolator interpolatorPhi(phiVector.size(), ROOT::Math::Interpolation::kAKIMA );
  interpolatorPhi.SetData(phiVector,phiValVector);
  
  int upsampleFactor=100;
  int interpCtr=0;
  double interpThetaArray[npointsTheta*upsampleFactor], interpThetaValArray[npointsTheta*upsampleFactor];
  double interpPhiArray[npointsPhi*upsampleFactor], interpPhiValArray[npointsPhi*upsampleFactor];
  
  for (double interp_theta=thetaArray[0];interp_theta<thetaArray[npointsTheta-1];
       interp_theta+=FINE_BIN_SIZE/double(upsampleFactor)){
    interpThetaArray[interpCtr]=interp_theta;
    interpThetaValArray[interpCtr]=interpolatorTheta.Eval(interp_theta);
    interpCtr++;
  }

  interpCtr=0;
  for (double interp_phi=phiArray[0];interp_phi<phiArray[npointsPhi-1];interp_phi+=FINE_BIN_SIZE/double(upsampleFactor)){
    interpPhiArray[interpCtr]=interp_phi;
    interpPhiValArray[interpCtr]=interpolatorPhi.Eval(interp_phi);
    interpCtr++;
  }


  //make a graph of
  TGraph *gInterpTheta=new TGraph((npointsTheta-1)*upsampleFactor,interpThetaArray,interpThetaValArray);
  TGraph *gInterpPhi=new TGraph((npointsPhi-1)*upsampleFactor,interpPhiArray,interpPhiValArray);

  /*
  if (drawFlag==1){ 
    TCanvas *cInterpolation=new TCanvas("cInterpolation","cInterpolation",800,800);
    cInterpolation->Divide(1,2);
    cInterpolation->cd(1);
    gInterpTheta->Draw("aP");
    gInterpTheta->SetMarkerStyle(22);
    gInterpTheta->SetMarkerSize(0.5);
    cInterpolation->cd(2);
    gInterpPhi->SetMarkerStyle(22);
    gInterpPhi->SetMarkerSize(0.5);
    gInterpPhi->Draw("aP");
//    cInterpolation->Print("interpolatedPeak.eps");

  }
*/

  //get the peak of the interpolated graph.
  Int_t peakBinTheta=FFTtools::getPeakBin(gInterpTheta);
  Double_t peakValTheta;
  Double_t peakTheta; 
  gInterpTheta->GetPoint(peakBinTheta,peakTheta,peakValTheta);
  Int_t peakBinPhi=FFTtools::getPeakBin(gInterpPhi);
  Double_t peakValPhi;
  Double_t peakPhi; 
  gInterpPhi->GetPoint(peakBinPhi,peakPhi,peakValPhi);

  delete gInterpTheta;
  delete gInterpPhi;
  
//  if (printFlag==1) cout<<"Peak of interpolation in theta: "<<peakTheta<<", and phi: "<<peakPhi<<endl;
  if (peakValTheta<peakValPhi) peakVal=peakValTheta;
  else peakVal=peakValPhi;

  //get the FWHM of the peak
  double halfPeakThetaLeft=0; 
  double halfPeakThetaRight=0;  
  
  for (int i=0;i<(npointsTheta-1)*upsampleFactor;i++){
    if (peakBinTheta-i>0){
      if (interpThetaValArray[peakBinTheta-i]<peakValTheta/2. 
          && interpThetaValArray[peakBinTheta-(i-1)]>peakValTheta/2.){
        halfPeakThetaLeft=(interpThetaArray[peakBinTheta-i]+interpThetaArray[peakBinTheta-(i-1)])/2.;
      }
    }
    if (halfPeakThetaLeft!=0) break;
  }
  //if (halfPeakThetaLeft==0) halfPeakThetaLeft=interpThetaArray[0];
  
  for (int i=0;i<(npointsTheta-1)*upsampleFactor;i++){
    if (peakBinTheta+i<(npointsTheta-1)*upsampleFactor){
      if (interpThetaValArray[peakBinTheta+i]<peakValTheta/2. 
          && interpThetaValArray[peakBinTheta+(i-1)]>peakValTheta/2.){
        halfPeakThetaRight=(interpThetaArray[peakBinTheta+i]+interpThetaArray[peakBinTheta+(i-1)])/2.;
      }    
    }
    if (halfPeakThetaRight!=0) break;
  }
  // if (halfPeakThetaRight==0) halfPeakThetaRight=interpThetaArray[(npointsTheta-1)*upsampleFactor];

  double halfPeakPhiLeft=0; 
  double halfPeakPhiRight=0;  
  for (int i=0;i<(npointsPhi-1)*upsampleFactor;i++){
    if (peakBinPhi-i>0){
      if (interpPhiValArray[peakBinPhi-i]<peakValPhi/2. && interpPhiValArray[peakBinPhi-(i-1)]>peakValPhi/2.){
        halfPeakPhiLeft=(interpPhiArray[peakBinPhi-i]+interpPhiArray[peakBinPhi-(i-1)])/2.;
      }
    }
    if (halfPeakPhiLeft!=0) break;
  }
  //if (halfPeakPhiLeft==0) halfPeakPhiLeft=interpPhiArray[0];
  
  for (int i=0;i<(npointsPhi-1)*upsampleFactor;i++){
    if (peakBinPhi+i<(npointsPhi-1)*upsampleFactor){
      if (interpPhiValArray[peakBinPhi+i]<peakValPhi/2. && interpPhiValArray[peakBinPhi+(i-1)]>peakValPhi/2.){
        halfPeakPhiRight=(interpPhiArray[peakBinPhi+i]+interpPhiArray[peakBinPhi+(i-1)])/2.;
      }    
    }
    if (halfPeakPhiRight!=0) break;
  }
  // if (halfPeakPhiRight==0) halfPeakPhiRight=interpPhiArray[(npointsPhi-1)*upsampleFactor];
  
  Double_t FWHMTheta, FWHMPhi; 
  if (halfPeakThetaRight==0 || halfPeakThetaLeft==0) FWHMTheta=-1;
  else FWHMTheta=halfPeakThetaRight-halfPeakThetaLeft;
  if (halfPeakPhiRight==0 || halfPeakPhiLeft==0) FWHMPhi=-1;
  else FWHMPhi=halfPeakPhiRight-halfPeakPhiLeft;
//  if (printFlag==1) cout<<"FWHM Theta: "<<FWHMTheta<<", FWHM Phi: "<<FWHMPhi<<endl;


  const double fwhm_conversion = 2 * sqrt(2 * log(2)); 

  peak->x = peakPhi; 
  peak->y = peakTheta; 
  peak->sigma_x = FWHMPhi / fwhm_conversion; 
  peak->sigma_y = FWHMTheta / fwhm_conversion; 
  peak->covar = 0;
  peak->val = peakVal; 
}



const int GAUS_A = 0; 
const int GAUS_X0 = 1; 
const int GAUS_Y0 = 2; 
const int GAUS_SIGMA_X = 3; 
const int GAUS_SIGMA_Y = 4; 
const int GAUS_COVAR = 5; 
const int GAUS_PEDESTAL = 6; 
const int GAUS_NPAR = 7; 

static double gaus2d(double *xx, double *ps)
{

  double A = ps[GAUS_A]; 
  double sx = ps[GAUS_SIGMA_X]; 
  double sy = ps[GAUS_SIGMA_Y]; 
  double x0 = ps[GAUS_X0]; 
  double y0 = ps[GAUS_Y0]; 
  double p =  ps[GAUS_COVAR]; 
  double z0 = ps[GAUS_PEDESTAL]; 

  double x= xx[0]; 
  double y= xx[1]; 

  return z0 + A / (2 * M_PI * sx * sy * sqrt(1-p*p)) 
         * exp( - (x-x0)*(x-x0) / (2 * sx * sx * (1-p*p)) 
                - (y-y0)*(y-y0) / (2 * sy * sy * (1-p*p))
                 + p * (x-x0) * (y-y0) / ( (1-p*p) * sx * sy)
              );
}

void UCorrelator::peakfinder::doPeakFindingGaussian(const TH2D * zoomed, FineMaximum * peak)
{


  TF2 gaus2dfn("fitter", gaus2d, GAUS_NPAR, 
                         zoomed->GetXaxis()->GetXmin(), zoomed->GetXaxis()->GetXmax(), 
                         zoomed->GetYaxis()->GetXmin(), zoomed->GetYaxis()->GetXmax() ); 

  int xmax, ymax, ignored; 
  int max_bin = zoomed->GetMaximumBin(xmax,ymax,ignored); 
  gaus2dfn.SetParameter(GAUS_A, zoomed->GetArray()[max_bin]); 
  gaus2dfn.SetParameter(GAUS_SIGMA_X, zoomed->GetRMS(1)); 
  gaus2dfn.SetParameter(GAUS_SIGMA_Y, zoomed->GetRMS(2)); 
  gaus2dfn.SetParameter(GAUS_X0, zoomed->GetXaxis()->GetBinCenter(xmax)); 
  gaus2dfn.SetParameter(GAUS_Y0, zoomed->GetYaxis()->GetBinCenter(ymax)); 
  gaus2dfn.SetParameter(GAUS_COVAR, 0); 
  gaus2dfn.SetParameter(GAUS_PEDESTAL, 0); 

  TH2D copy = *zoomed;
  copy.Fit(&gaus2dfn, "NQ"); 

  peak->sigma_x = gaus2dfn.GetParameter(GAUS_SIGMA_X); 
  peak->sigma_y = gaus2dfn.GetParameter(GAUS_SIGMA_Y); 
  peak->covar = gaus2dfn.GetParameter(GAUS_COVAR) * peak->sigma_x * peak->sigma_y; 
  peak->x = gaus2dfn.GetParameter(GAUS_X0); 
  peak->y = gaus2dfn.GetParameter(GAUS_Y0); 
  peak->val = gaus2dfn.GetParameter(GAUS_A) / (2 * M_PI * peak->sigma_x * peak->sigma_y * sqrt(1-peak->covar *peak->covar)); 
}


void UCorrelator::peakfinder::doInterpolationPeakFindingBicubic(const TH2D * zoomed, FineMaximum * peak)
{
  fprintf(stderr,"doInterpolationPeakFindingBicubic not implemented yet!\n"); 


}



