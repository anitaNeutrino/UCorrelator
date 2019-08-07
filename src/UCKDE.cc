#include "UCKDE.h" 
#include "TMath.h" 

ClassImp(UCorrelator::KDE2D); 


static double dist2(double x0, double y0, double x1, double y1) 
{

  return (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1); 

}

//naive algorithm for now... probably very slow! 
static int nearest_neighbors(int n, double x, double y, int Npts, const double * vx, 
                            const double * vy, int * nearest_i, double * nearest_d2) 
{
  
  int nfound = 0; 
  
  for (int i = 0; i < Npts; i++) 
  {

    double d2 = dist2(x,y,vx[i],vy[i]); 

    if (nfound == 0) 
    {
      nearest_d2[0]=d2; 
      nearest_i[0] = i; 
      nfound++; 
    }
    else if (d2 < nearest_d2[nfound-1]) 
    {

      int last = nfound-1;
      //if we are below capacity, increase capacity! 
      if (nfound < n) nfound++; 
      
      int new_index = last; 
      for (int ii = 0 ; ii < last-1; ii++) 
      {
        if (d2 < nearest_d2[ii]) 
        {
          new_index = ii; 
          break; 
        }
      }

      //now shift things over 1

      for (int ii = nfound-1; ii > new_index; ii--) 
      {
        nearest_d2[ii] = nearest_d2[ii-1]; 
        nearest_i[ii] = nearest_i[ii-1]; 
      }
      
      nearest_d2[new_index] = d2; 
      nearest_i[new_index] = i; 
    }
  }
  return nfound; 
}



UCorrelator::KDE2D::KDE2D(int n, const double *xin, const double * yin, const double * weights, KDE2DOptions opt) 
  : x(xin, xin+n), y(yin,yin+n) 
{
  if (weights) 
  {
    w.insert(w.end(), weights, weights+n); 
  }

  double sx = opt.sigma_x ? opt.sigma_x: pow(n,-1./6)*TMath::RMS(n,xin);
  double sy = opt.sigma_y ? opt.sigma_y: pow(n,-1./6)*TMath::RMS(n,yin);

  // sigma_x/y will be of size 1 if not adaptive. 
  
  if (opt.adaptiveness ==  0) 
  {

    sigma_x.push_back(sx); 
    sigma_y.push_back(sy); 
    rho.push_back(opt.corr); 
    double sum_of_weights = n; 
    if (w.size() >  0) 
    {
      sum_of_weights = 0; 
      for (unsigned i = 0; i < w.size(); i++) sum_of_weights += w[i]; 
    }
    W = 1./( 2*sum_of_weights*M_PI * sx*sy*sqrt(1-opt.corr*opt.corr));   //check this
  }

  else 
  {
    sigma_x.reserve(n);
    sigma_y.reserve(n);
    rho.reserve(n);

    std::vector<double> nearest_d2(opt.adaptiveness+1); 
    std::vector<int> nearest_i(opt.adaptiveness+1); 
    double Winv =  0; 
    for (int i = 0; i < n; i++) 
    {
      int found = nearest_neighbors(opt.adaptiveness+1, xin[i],yin[i], n, xin, yin, &nearest_i[0], &nearest_d2[0]); 

      double d2 =nearest_d2[found-1]; 
      double d = sqrt(d2); 
      sigma_x.push_back(sx * d); 
      sigma_y.push_back(sy * d); 
      rho.push_back(opt.corr * d); 
      double this_w = w.size() ? w[i] : 1; 
      Winv += this_w/( 2*d2*M_PI * sx*sy*sqrt(1-d2*opt.corr*opt.corr));  //check this 
    }
    W = 1./Winv; 
  } 
}


//TODO vectorize this 
double UCorrelator::KDE2D::operator()(double xt, double yt) const
{

  double sxinv = 1./sigma_x[0];
  double syinv = 1./sigma_y[0];
  double sxinv2 = sxinv*sxinv;
  double syinv2 = syinv*syinv; 
  double r = rho[0];
  double one_minus_r_inv = 1./(1-r);
  double sxyinv = r == 0 ? 0 : sxinv*syinv; 

  double sum = 0;
  for (unsigned i = 0; i < x.size(); i++) 
  {
    if (sigma_x.size() > 1 && i > 0) 
    {
      sxinv = 1./sigma_x[i]; 
      syinv = 1./sigma_y[i]; 
      sxinv2 = sxinv*sxinv;
      syinv2 = syinv*syinv; 
      r = rho[i]; 
      one_minus_r_inv = 1./(1-r);
      sxyinv = r == 0 ? 0 : sxinv*syinv; 
    }

    double xd = ( x[i] - xt);
    double yd = ( y[i] - yt);
    double z = xd*xd*sxinv2 + yd*yd*syinv2; 
    if (r!=0) z -= 2*r*xd*yd*sxyinv; 
    double this_w = w.size() == 0 ? 1 : w[i]; 
    sum += this_w * (r == 0 ?  exp ( - 0.5* z ) : exp(-0.5*z*one_minus_r_inv)); 
  }

  return sum *W; 
}


void UCorrelator::KDE2D::getNSigmaBounds(double nsig, double & xmin, double & xmax, double & ymin, double & ymax) const
{

  double sig_x = sigma_x[0]; 
  double sig_y = sigma_y[0]; 

  xmin = x[0] - sig_x * nsig; 
  xmax = x[0] + sig_x * nsig; 
  ymin = y[0] - sig_y * nsig; 
  ymax = y[0] + sig_y * nsig; 

  for (unsigned i = 1; i <  x.size() ; i++) 
  {
    if (sigma_x.size() > 1) 
    {
      sig_x = sigma_x[i]; 
      sig_y = sigma_y[i]; 
    }

    double test_xmin = x[i] - sig_x * nsig; 
    double test_xmax = x[i] + sig_x * nsig; 
    double test_ymin = y[i] - sig_y * nsig; 
    double test_ymax = y[i] + sig_y * nsig; 


    if (test_xmin < xmin) xmin = test_xmin; 
    if (test_ymin < ymin) ymin = test_ymin; 
    if (test_xmax > xmax) xmax = test_xmax; 
    if (test_ymax > ymax) ymax = test_ymax; 
  }
}


TF2 * UCorrelator::KDE2D::makeTF2(const char * name, double nsigma_bounds) const
{

  double xmin,xmax,ymin,ymax; 
  getNSigmaBounds(nsigma_bounds,xmin,xmax,ymin,ymax); 

  return new TF2(name, this, xmin,xmax,ymin,ymax, 0,2);
}

TH2 * UCorrelator::KDE2D::makeTH2(double binwidth_x, double binwidth_y, const char * name, const char * title, double nsigma_bounds, double auto_factor) const
{

  double xmin,xmax,ymin,ymax;
  getNSigmaBounds(nsigma_bounds,xmin,xmax,ymin,ymax); 

  if (!binwidth_x) 
  {
    binwidth_x = sigma_x[0] * auto_factor; 
    for (unsigned i = 1; i < sigma_x.size(); i++) binwidth_x = TMath::Min(binwidth_x, sigma_x[i]*auto_factor);
  }

  if (!binwidth_y) 
  {
    binwidth_y = sigma_y[0] * auto_factor; 
    for (unsigned i = 1; i < sigma_y.size(); i++) binwidth_y = TMath::Min(binwidth_y, sigma_y[i]*auto_factor);
  }

  int nbinsx = ceil((xmax-xmin)/binwidth_x); 
  int nbinsy = ceil((ymax-ymin)/binwidth_x); 
  xmax = xmin + nbinsx * binwidth_x; 
  ymax = ymin + nbinsy * binwidth_y; 

  TH2 * h = new TH2D ( name,title, nbinsx, xmin, xmax, nbinsy, ymin, ymax); 

  h->SetEntries(x.size()); 
  for (int i = 1; i <= h->GetNbinsX(); i++) 
  {
    for (int j = 1;  j <= h->GetNbinsY(); j++) 
    {
      h->SetBinContent(i,j, this->operator()(h->GetXaxis()->GetBinCenter(i), h->GetYaxis()->GetBinCenter(j))); 
    }
  }
  return h; 
}
