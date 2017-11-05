#include "PeakFinder.h"
#include "TH2.h" 
#include "TF2.h" 
#include "FFTtools.h"
#include "Minuit2/Minuit2Minimizer.h" 
#include "Math/Interpolator.h"
#include "TGraph.h"

#ifdef UCORRELATOR_USE_EIGEN_FOR_PEAK_FINDER
#include <Eigen/Dense> 
#else 
#include "TMatrix.h"
#include "TDecompSVD.h"
#endif 

#include "TError.h" 

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

static int findMaximum(const TH2D* hist, std::vector<char> * notallowed)
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

static void maskNearbyBins(const TH2D * hist, double distance, int bin, std::vector<char> * used)
{

  const int width = hist->GetNbinsX()+2; 
  const int height = hist->GetNbinsY()+2; 
  int int_dist = int(distance+0.5); 

  double dw = hist->GetXaxis()->GetBinWidth(1); 
  double dh = hist->GetYaxis()->GetBinWidth(1); 
  double dw2 = dw*dw; 
  double dh2 = dh*dh; 


  int i0 = bin % width;
  int j0 = bin / width; 

  double dist2 = distance*distance; 

  for (int jj = -int_dist; jj <= int_dist; jj++)
  {
    int j = j0 +jj; 
    if (j < 1) continue; 
    if (j >= height) continue; 
    int j2 = jj*jj; 
    for (int ii = -int_dist; ii <= int_dist; ii++) 
    {
      int i = i0 +ii; 
      if (i < 1) i += width; 
      if (i >= width) i-=width; 
      if (ii*ii*dw2 + j2*dh2 > dist2) continue; 

      int bin2use =  i+ j * width; 
//      printf("Disallowing %d\n", bin2use); 
      (*used)[bin2use] = 1; 
    }
  }
}


int UCorrelator::peakfinder::findIsolatedMaxima(const TH2D * hist, double distance, int Nmaxima, RoughMaximum * maxima, double minPhi, double maxPhi, double minTheta, double maxTheta, bool exclude, bool use_bin_center)
{
  int width = hist->GetNbinsX()+2; 
  int height = hist->GetNbinsY()+2; 
  std::vector<char> used( width*height, false); 

  int nfound = 0; 
	
	std::vector<int> row_not_allowed(0);
	std::vector<int> col_not_allowed(0);
	double minWrap = 360.;
	double maxWrap = 0.;
	if(minPhi < 0) minWrap += minPhi;
	if(maxPhi > 360) maxWrap -= 360;

	//if min and maxes are set, blocks out corresponding rows/columns
	if(minTheta || maxTheta)
	{
		for(int ybin = 2; ybin < hist->GetNbinsY(); ybin++)
		{
			double yCenter = hist->GetYaxis()->GetBinCenter(ybin);
			if(yCenter <= minTheta || yCenter >= maxTheta) row_not_allowed.push_back(ybin);
		}
	}
	if(minPhi || maxPhi)
	{
		for(int xbin = 1; xbin < hist->GetNbinsX()+1; xbin++)
		{
			double xCenter = hist->GetXaxis()->GetBinCenter(xbin);
			if((xCenter <= minPhi && xCenter >= maxWrap) || (xCenter >= maxPhi && xCenter <= minWrap))
			{
				col_not_allowed.push_back(xbin);
			}
		}
	}

  for (unsigned i = 0; i < row_not_allowed.size(); i++)
  {
    memset(&used[(width)*row_not_allowed[i]], 1, width); 
  }
	for (unsigned i = 0; i < col_not_allowed.size(); i++)
  {
		for(unsigned j = 0; j < (unsigned) height; j++) used[col_not_allowed[i] + (j*width)] = 1;
  }

	//if set to exclude, this flips all 1s and 0s (max sure phi and theta are actually set if using this or it will just block out everything)
	if(exclude)
	{
		for(unsigned i = 0; i <  (unsigned)width*height; i++) used[i] = used[i] xor 1;
	}
  
	//block out bins closest to top and bottom since we don't want maxima on the top or bottom edge 

  int rows_not_allowed[] = {1,hist->GetNbinsY()}; 
	for(unsigned i = 0; i < sizeof(rows_not_allowed) / sizeof(*rows_not_allowed); i++)
	{
		memset(&used[(width)*rows_not_allowed[i]], 1, width);
	}

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
  p->chisq = chisq; 
}

const int X2COEFF = 0; 
const int Y2COEFF = 1; 
const int XYCOEFF = 2; 
const int XCOEFF = 3; 
const int YCOEFF = 4; 
const int CCOEFF = 5; 



#ifdef UCORRELATOR_USE_EIGEN_FOR_PEAK_FINDER

template <unsigned N> 
static const Eigen::Matrix<double,N*N,6> & getM()
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
  return M; 
}

template <unsigned N> 
static const Eigen::JacobiSVD<Eigen::Matrix<double, N*N,6> > & getSVD()
{
  static Eigen::JacobiSVD<Eigen::Matrix<double, N*N,6> > svd(getM<N>(),Eigen::ComputeFullU | Eigen::ComputeFullV); 
  return svd; 
}
#else

template <unsigned N>
static const TMatrixD & getM()
{
  int halfway = N/2; 
  static TMatrixD M(N*N,6); 
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

  return M;
}

template <unsigned N> 
static const TDecompSVD & getSVD()
{
  static TDecompSVD svd(getM<N>()); 
  svd.Decompose(); 
  return svd; 


}

#endif


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
#ifdef UCORRELATOR_USE_EIGEN_FOR_PEAK_FINDER
  static const Eigen::JacobiSVD<Eigen::Matrix<double, N*N, 6> > & svd = getSVD<N>();  //I  hope this works! 
  Eigen::Matrix<double,N*N,1> Z; 
#else
  static const TDecompSVD  & svd = getSVD<N>(); 
  TVectorD Z(N*N); 

#endif

  // double xcenter = hist->GetXaxis()->GetBinCenter(xmax); 
  // double ycenter = hist->GetYaxis()->GetBinCenter(ymax); 
  //peng, change the center of bin to low edge of bin.
  double xcenter = hist->GetXaxis()->GetBinLowEdge(xmax);
  double ycenter = hist->GetYaxis()->GetBinLowEdge(ymax); 

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
#ifdef UCORRELATOR_USE_EIGEN_FOR_PEAK_FINDER
  Eigen::Matrix<double,6,1> B = svd.solve(Z); 
#else
  bool ok = true; 
  TVectorD B = ((TDecompSVD & )svd).Solve(Z,ok);  // I think this is effectively const since the matrix is fixed, but maybe it'll backfire with multiple threads. In that case, use EIGEN I guess?  
  if (!ok) 
  {
    fprintf(stderr, "Warning in doQuadraticPeakFinding: SVD decomposition failed!\n"); 
  }
#endif
    
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
  peak->sigma_x =  sqrt(-2*b*peak->val / fabs(4*a*b - c*c))*dx; //i don't think that fabs is correct
  peak->sigma_y =  sqrt(-2*a*peak->val / fabs(4*a*b - c*c))*dy;
  peak->covar = c*peak->val / fabs(4*a*b - c*c)*dx*dy; 
#ifdef UCORRELATOR_USE_EIGEN_FOR_PEAK_FINDER
  peak->chisq = (getM<N>() * B - Z).squaredNorm()/ (N*N); 
#else
  peak->chisq = (getM<N>() * B - Z).Norm2Sqr()/ (N*N); 
#endif

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
  peak->chisq = 0; // not sure what the best way to do this is 

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
  peak->chisq = gaus2dfn.GetChisquare()/ gaus2dfn.GetNDF(); 
}


const double  M1_elems[] = { 1,0,0,0,
                             0,0,1,0,
                             -3,3,-2,-1,
                             2,-2,1,1}; 

const static TMatrixD M1(4,4, M1_elems); 
const double  M2_elems[] = { 1,0,-3,2,
                             0,0,3,-2,
                             0,1,-2,1,
                             0,0,-1,1}; 

const static TMatrixD M2(4,4, M2_elems); 




class NegativeBicubicFunction : public ROOT::Math::IGradientFunctionMultiDim
{

  public:
    NegativeBicubicFunction( const double m[4][4], double xwidth, double ywidth) ; 
    NegativeBicubicFunction( const NegativeBicubicFunction & other) { memcpy(a,other.a, sizeof(a)); } 

    virtual double DoEval(const double * x) const; 
    virtual double DoDerivative(const double * x, unsigned int coord) const; 
    unsigned NDim() const { return 2; } 

    virtual ROOT::Math::IBaseFunctionMultiDim * Clone() const { return new NegativeBicubicFunction(*this); }

  private: 
    double a[4][4]; 
}; 


NegativeBicubicFunction::NegativeBicubicFunction (const double m[4][4], double xw, double yw) 
{
  TMatrixD Mf(4,4); 
  //construct the bicubic polynomial over this bin 
      
  // This uses notation similar to wikipedia  


  double f_00 = m[1][1]; 
  double f_11 = m[2][2]; 
  double f_01 = m[1][2]; 
  double f_10 = m[2][1]; 


  double fx_00 = (m[2][1] - m[0][1]) / (2*xw); 
  double fx_10 = (m[3][1] - m[1][1]) / (2*xw); 
  double fx_01 = (m[2][2] - m[0][2]) / (2*xw); 
  double fx_11 = (m[3][2] - m[1][2]) / (2*xw); 
      
  double fy_00 = (m[1][2] - m[1][0]) / (2*yw); 
  double fy_10 = (m[2][2] - m[2][0]) / (2*yw); 
  double fy_01 = (m[1][3] - m[1][1]) / (2*yw); 
  double fy_11 = (m[2][3] - m[2][1]) / (2*yw); 
 
  double fxy_00 = ( (m[2][2] - m[0][2]) -(m[2][0] - m[0][0])) / (4*xw*yw) ; 
  double fxy_01 = ( (m[2][3] - m[0][3]) -(m[2][1] - m[0][1])) / (4*xw*yw) ; 
  double fxy_10 = ( (m[3][2] - m[1][2]) -(m[3][0] - m[1][0])) / (4*xw*yw) ; 
  double fxy_11 = ( (m[3][3] - m[1][3]) -(m[3][1] - m[1][1])) / (4*xw*yw) ; 

  Mf(0,0) = f_00; Mf(0,1) = f_01; Mf(0,2) = fy_00; Mf(0,3) = fy_01; 
  Mf(1,0) = f_10; Mf(1,1) = f_11; Mf(1,2) = fy_10; Mf(1,3) = fy_11; 
  Mf(2,0) = fx_00; Mf(2,1) = fx_01; Mf(2,2) = fxy_00; Mf(2,3) = fxy_01; 
  Mf(3,0) = fx_10; Mf(3,1) = fx_11; Mf(3,2) = fxy_10; Mf(3,3) = fxy_11; 

 //now there is probably some smart way to reuse some of the information between bins, but I dont' feel like figuring it out right now 

  TMatrix Ma = M1 * Mf * M2; 
  for (int i = 0; i < 4; i++) 
  {
    for (int j = 0; j < 4; j++)
    {
      a[i][j] = Ma(i,j); 
    }
  }

}

double NegativeBicubicFunction::DoEval(const double * p) const 
{
  double x[4]; 
  double y[4]; 
  x[0] = 1; 
  y[0] = 1; 
  for (int i = 1; i < 3; i++) 
  {
    x[i] = p[0] * x[i-1]; 
    y[i] = p[1] * y[i-1]; 
  }


  double val = 0; 

  for (int i = 0; i < 3; i++) 
  {
    for (int j = 0; j < 3; j++) 
    {
      val += a[i][j] * x[i] * y[j]; 
    }
  }

  return -val; 
}

double NegativeBicubicFunction::DoDerivative(const double * p, unsigned int coord) const
{
  double x[4]; 
  double y[4]; 
  x[0] = 1; 
  y[0] = 1; 
  for (int i = 1; i < 3; i++) 
  {
    x[i] = p[0] * x[i-1]; 
    y[i] = p[1] * y[i-1]; 
  }


  double val = 0; 


  if (coord == 0) 
  {
    for (int i = 1; i < 3; i++) 
    {
      for (int j = 0; j < 3; j++) 
      {
        val += a[i][j] * x[i-1] * y[j] * i; 
      }
    }
  }
  else 
  {
    for (int i = 0; i < 3; i++) 
    {
      for (int j = 1; j < 3; j++) 
      {
        val += a[i][j] * x[i] * y[j-1]*j; 
      }
    }
  }

  return -val; 
}


static __thread ROOT::Minuit2::Minuit2Minimizer *min = 0; 

void UCorrelator::peakfinder::doInterpolationPeakFindingBicubic(const TH2D * zoomed, FineMaximum * peak, double min_to_consider )
{

  TMatrixD Mf(4,4); 
  double xwidth = zoomed->GetXaxis()->GetBinWidth(1); 
  double ywidth = zoomed->GetYaxis()->GetBinWidth(1); 
  double m[4][4]; 

  if (!min) min = new ROOT::Minuit2::Minuit2Minimizer; 

  double stupid_max = zoomed->GetMaximum(); 

  peak->val = 0; 

  for (int i = 2; i < zoomed->GetNbinsX(); i++) 
  {
    for (int j = 2; j < zoomed->GetNbinsY(); j++) 
    {
      //stupid heuristic quick fail 
      bool ok = false; 
      for (int ii = -1; ii < 3; ii++)
      {
        for (int jj = -1; jj < 3; jj++)
        {
          double val = zoomed->GetBinContent(i+ii, j+ii); 
          if (val > stupid_max * min_to_consider) ok = true; 

          m[ii+1][jj+1] = zoomed->GetBinContent(i+ii,j+jj); 
        }
      }

      if (!ok) continue; 

      NegativeBicubicFunction f(m,xwidth,ywidth); 

      min->Clear(); 
      min->SetFunction(f); 
      min->SetLimitedVariable(0,"x", 0.5, 0.1, 0,1); 
      min->SetLimitedVariable(1,"y", 0.5, 0.1, 0,1); 


      //shut it up 
      int old_level = gErrorIgnoreLevel; 
      min->Minimize(); 
      gErrorIgnoreLevel = old_level; 


      if (-min->MinValue() > peak->val)
      {
        peak->val = -min->MinValue(); 
        peak->x = min->X()[0]; 
        peak->y = min->X()[1]; 
        peak->sigma_x = sqrt(min->CovMatrix(0,0)); 
        peak->sigma_y = sqrt(min->CovMatrix(1,1)); 
        peak->covar = min->CovMatrix(0,1); 
        peak->chisq = 0; 
      }
    }
  }
}


void UCorrelator::peakfinder::doPeakFindingHistogram(const TH2D* hist, FineMaximum * peak) 
{

  int max_bin = hist->GetMaximumBin(); 
  int bx,by,bz; 
  hist->GetBinXYZ(max_bin,bx,by,bz); 

  peak->val = hist->GetBinContent(bx,by); 
  peak->x =  hist->GetXaxis()->GetBinCenter(bx); 
  peak->y =  hist->GetYaxis()->GetBinCenter(by); 

  //compute covariance matrix around point 
  
  double var_x = 0;
  double var_y = 0; 
  double covar = 0; 
  double total = 0; 

  for (int i = 1; i <= hist->GetNbinsX(); i++)
  {
    for (int j = 1; j <= hist->GetNbinsY(); j++)
    {
      double val = hist->GetBinContent(i,j); 
      double x = hist->GetXaxis()->GetBinCenter(i); 
      double y = hist->GetYaxis()->GetBinCenter(j); 

      total += val; 
      var_x += val  * (x-peak->x) * (x-peak->x); 
      var_y += val  * (y-peak->y) * (y-peak->y); 
      covar += val  * (y-peak->y) * (x-peak->x); 
    }
  }


  peak->sigma_x = sqrt( var_x / (total-1)); 
  peak->sigma_y = sqrt( var_y / (total-1)); 
  peak->covar =  covar / (total - 1); 
  peak->chisq = 0;  
}

