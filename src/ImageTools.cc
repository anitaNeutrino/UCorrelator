#include "TProfile2D.h"
#include "TH3.h" 
#include "UCImageTools.h" 
#include "TString.h" 
#include <algorithm>
#include <set>
#include <math.h>

static TH1* makeTH1(char type, const char * name, const char * title, int nbins, double xmin, double xmax) 
{
  switch(type)
  {
    case 'C': 
      return new TH1C(name,title,nbins,xmin,xmax); 
    case 'S': 
      return new TH1S(name,title,nbins,xmin,xmax); 
    case 'I': 
      return new TH1I(name,title,nbins,xmin,xmax); 
    case 'F': 
      return new TH1F(name,title,nbins,xmin,xmax); 
    case 'D': 
    default:
      return new TH1D(name,title,nbins,xmin,xmax); 
  }
}






TProfile2D *UCorrelator::TH2toTProfile2D(TH2* inTH2) {

  int nBinX = inTH2->GetNbinsX();
  int nBinY = inTH2->GetNbinsY();
  int xMin  = inTH2->GetXaxis()->GetBinLowEdge(1);
  int yMin  = inTH2->GetYaxis()->GetBinLowEdge(1);
  int xMax  = inTH2->GetXaxis()->GetBinUpEdge(nBinX);
  int yMax  = inTH2->GetYaxis()->GetBinUpEdge(nBinY);
  TProfile2D *mapProfile = new TProfile2D("temp","temp",nBinX,xMin,xMax,nBinY,yMin,yMax);

  return mapProfile;

}

void UCorrelator::fillTProfile2DWithTH2(TProfile2D *prof, TH2* hist) {

  int nBinX = hist->GetNbinsX();
  int nBinY = hist->GetNbinsY();

  for (int binX=0; binX<nBinX+1; binX++) {
    double x = hist->GetXaxis()->GetBinCenter(binX+1);
    for (int binY=0; binY<nBinY+1; binY++) {
      double y = hist->GetYaxis()->GetBinCenter(binY+1);
          prof->Fill(x,y,hist->GetBinContent(binX,binY));
    }
  }

  return;
}



TH2* UCorrelator::rotateHistogram(const TH2* inHist, double rotate) {
  
  //get all the info about the graph
  int nBinX = inHist->GetNbinsX();
  int nBinY = inHist->GetNbinsY();
  double xMin  = inHist->GetXaxis()->GetBinLowEdge(1);
  double yMin  = inHist->GetYaxis()->GetBinLowEdge(1);
  double xMax  = inHist->GetXaxis()->GetBinUpEdge(nBinX);
  double yMax  = inHist->GetYaxis()->GetBinUpEdge(nBinY);

  const char* title = inHist->GetTitle();

  TH2D *outHist = new TH2D("rotateHistogram",title,nBinX,xMin,xMax,nBinY,yMin,yMax);

  for (Int_t iX = 0; iX < nBinX; iX++) {

    double rotX = inHist->GetXaxis()->GetBinCenter(iX+1) - rotate;
    //get that rotated value into the range
    while (rotX < 0) rotX += 360.;
    while (rotX >= 360) rotX -= 360.;

    for (Int_t iY = 0; iY < nBinY; iY++) {     

      double Y = inHist->GetYaxis()->GetBinCenter(iY+1);

      double value = inHist->GetBinContent(iX+1,iY+1);

      outHist->Fill(rotX,Y,value);
    }
  }

  return outHist;
}




double UCorrelator::getZRMS(const TH2* hist, const int * lims) {

  int numBinsX = hist -> GetNbinsX();
  int numBinsY = hist -> GetNbinsY();

  int minBinY = lims ? lims[2] : 1; 
  int maxBinY = lims ? lims[3] : numBinsY; 

  double mean = 0, sqMean = 0;
  int numBinsUsedX = 0; 

  for (int binX = 1; binX <= numBinsX; ++binX) {  //  For histograms, index 0 and length + 1 represent underflow and overflow bins, respectively.
    if (lims)//slow more complicated code for now
    {
      if (lims[0] < lims[1] && (binX < lims[0] || binX > lims[1])) continue;
      if (lims[0] > lims[1] && (binX > lims[0] || binX < lims[1])) continue;
      numBinsUsedX++; 
    }
    for (int binY = minBinY; binY <= maxBinY; ++binY) {

      double binValue = hist -> GetBinContent(binX, binY);
      mean += binValue;
      sqMean += binValue * binValue;
    }
  }
  int numBinsUsedY = lims ? (maxBinY-minBinY + 1) : numBinsY; 
  
  if (!lims) numBinsUsedX = numBinsX;

  mean /= numBinsUsedX * numBinsUsedY;
  sqMean /= numBinsUsedX * numBinsUsedY;

  return sqrt(sqMean - mean * mean);
}



static const TAxis * getAxis(const TH2* H, int axis) 
{
  switch(axis) 
  {

    case 1: 
      return H->GetXaxis(); 
    case 2: 
      return H->GetYaxis(); 
    default:
      fprintf(stderr,"UCorrelator::image: axis must be 1 (x) or 2 (y)\n"); 
      return 0; 
  }
}


TH1* UCorrelator::image::getPctileProjection(const TH2 * H, int axis, double pct, bool ignoreEmpty, const TH1 * norm) 
{
  TString name; 
  TString title; 
  name.Form("%s_%g_proj_%d", H->GetName(), pct, axis); 
  title.Form("%s %g%% projection of axis %d", H->GetTitle(), 100*pct, axis); 
  
  TH1 * h = makeTH1(H->ClassName()[3], name.Data(), title.Data(), getAxis(H,axis)->GetNbins(), getAxis(H,axis)->GetXmin(), getAxis(H,axis)->GetXmax()); 

  int nOrth = getAxis(H,3-axis)->GetNbins(); 

  std::set<int> ignore; 

  if (ignoreEmpty) 
  {
    for (int j = 1; j <= nOrth; j++) 
    {
      if (norm && norm->GetBinContent(j) == 0) 
      {
        ignore.insert(j); 
      }
      else
      {
        bool empty = true; 

        for (int i =1; i <= h->GetNbinsX(); i++) 
        {
          if( H->GetBinContent( axis == 1 ? i : j, axis == 1 ? j : i))
          {
            empty = false; 
            break; 
          }
        }

        if (empty) 
          ignore.insert(j); 
        }
    }

  }

  for (int i = 1; i <= h->GetNbinsX(); i++) 
  {
    std::vector<double> v;
    v.reserve(nOrth); 
    for (int j = 1; j <= nOrth; j++)
    {
      if (ignore.count(j)) 
        continue; 

      v.push_back( H->GetBinContent( axis == 1 ? i : j, axis == 1 ? j : i)); 
    }
    std::nth_element(v.begin(), v.begin() + (v.size()-1) * pct, v.end()); 
    h->SetBinContent(i, v[(v.size()-1)*pct]); 
  }

  return h; 
} 


static double bicubicFunction(double x, double * params)
{
  double * a = params+1;  //hack to use same notation as in wikipedia :)
  return 0.5 * (2*a[0] + x*(-a[-1] + a[1]) + x*x*(2*a[-1]-5*a[0]+4*a[1] -a[2]) + x*x*x*(-a[-1] + 3*a[0] - 3*a[1] + a[2])); 
}


double UCorrelator::image::interpolate(const TH2 *h, double x, double y, InterpolationType type, InterpolationEdgeBehavior flags_x, InterpolationEdgeBehavior flags_y, bool centers)
{

  if (type == NEAREST) 
  {
    int xbin = h->GetXaxis()->FindFixBin(x); 
    int ybin = h->GetXaxis()->FindFixBin(x); 

    if (centers)  //this is easy 
    {
      return h->GetBinContent(xbin,ybin); 
    }

    //otherwise, we want to find the closest corner. 

    if (x > h->GetXaxis()->GetBinCenter(xbin))
    {
      xbin = xbin < h->GetNbinsX() || flags_x == EXTEND ? xbin+ 1  :  flags_x == PERIODIC ? 1  : xbin; 
    }

    if (y > h->GetYaxis()->GetBinCenter(ybin))
    {
      ybin = ybin < h->GetNbinsY() || flags_y == EXTEND ? ybin+ 1  :  flags_y == PERIODIC ? 1  : ybin; 
    }

    return h->GetBinContent(xbin,ybin); 
  }


  if (type == BILINEAR)  
  {

    //find the bin
    
    int xbin = h->GetXaxis()->FindFixBin(x); 
    int ybin = h->GetYaxis()->FindFixBin(y); 

    if (centers && x < h->GetXaxis()->GetBinCenter(xbin) )
      xbin = (xbin == 1 && flags_x == PERIODIC) ? h->GetNbinsX() : xbin-1; 

    if (centers && y < h->GetYaxis()->GetBinCenter(ybin) )
      ybin = (ybin == 1 && flags_y == PERIODIC) ? h->GetNbinsY() : ybin-1; 



    //bins used to get the value 
    int xbin2 = xbin < h->GetNbinsX() ? xbin + 1 : 
                flags_x == PERIODIC   ? 1        :
                flags_x == EXTEND     ? xbin     :
                xbin+1; //this will be the overflow bin which I suppose will usually be zero 

    int ybin2 = ybin < h->GetNbinsY() ? ybin + 1 : 
                flags_y == PERIODIC   ? 1        :
                flags_y == EXTEND     ? ybin     :
                ybin+1; //this will be the overflow bin which I suppose will usually be zero 

    //get the values 
    double q11  = h->GetBinContent(xbin,ybin); 
    double q21  = h->GetBinContent(xbin2,ybin); 
    double q12  = h->GetBinContent(xbin,ybin2); 
    double q22  = h->GetBinContent(xbin2,ybin2); 

    double xwidth = centers ? h->GetXaxis()->GetBinCenter(xbin+1) - h->GetXaxis()->GetBinCenter(xbin) : h->GetXaxis()->GetBinWidth(xbin); 
    double ywidth = centers ? h->GetYaxis()->GetBinCenter(ybin+1) - h->GetYaxis()->GetBinCenter(ybin) : h->GetYaxis()->GetBinWidth(ybin); 

    // TODO: special case the uniform bin size case
    
    double xrel = x - (centers ? h->GetXaxis()->GetBinCenter(xbin) :  h->GetXaxis()->GetBinLowEdge(xbin)); 
    double yrel = y - (centers ? h->GetYaxis()->GetBinCenter(ybin) :  h->GetYaxis()->GetBinLowEdge(ybin)); 
    double xrel2 = xwidth - xrel; 
    double yrel2 = ywidth - yrel; 

    return 1./ (xwidth * ywidth)  * ( q11 * xrel2 * yrel2 + q21 *xrel * yrel2 + q12 * xrel2 * yrel  + q22 * xrel * yrel); 
  }

  if (type  == BICUBIC)
  {

    if (h->GetXaxis()->GetXbins()->fN || h->GetYaxis()->GetXbins()->fN)
    {
      fprintf(stderr,"WARNING, bicubic interpolation does not work for non-uniform binning. Reverting to bilinear.\n"); 
      return interpolate(h,x,y,BILINEAR,flags_x,flags_y); 
    }

    
    double xmin = h->GetXaxis()->GetXmin(); 
    double xmax = h->GetXaxis()->GetXmax(); 
    double ymin = h->GetYaxis()->GetXmin(); 
    double ymax = h->GetYaxis()->GetXmax(); 
    double xwidth = (xmax-xmin) / (h->GetNbinsX()); 
    double ywidth = (ymax-ymin) / (h->GetNbinsY()); 

  
    int binx0 = (int) ((x - xmin) / xwidth + centers ? 0.5 : 0); 
    int biny0 = (int) ((y - ymin) / ywidth + centers ? 0.5 : 0 ); 


    double b[4]; 
    for (int iy = -1; iy <3; iy++)
    {
      double a[4]; 
      int biny = biny0 + iy; 

      if (biny < 1) 
      { 
        biny = flags_y == PERIODIC ? h->GetNbinsY() + biny : 
               flags_y == EXTEND  ? 1 : 
               0 ;
      }
      
      if (biny > h->GetNbinsY()  ) 
      {

        biny = flags_y == PERIODIC ? biny - h->GetNbinsY() : 
               flags_y == EXTEND ? h->GetNbinsY()          : 
               h->GetNbinsY() + 1;

      }

      for (int ix = -1; ix <3; ix++)
      {
        int binx = binx0 + ix; 

        if (binx < 1) 
        { 
          binx = flags_x == PERIODIC ? h->GetNbinsX() + binx : 
                 flags_x == EXTEND  ? 1 : 
                 0 ;
        }
        
        if (binx> h->GetNbinsX()  ) 
        {

          binx = flags_x == PERIODIC ? binx - h->GetNbinsX() : 
                 flags_x == EXTEND ? h->GetNbinsX()          : 
                 h->GetNbinsX() + 1;

        }


        a[ix+1] = h->GetBinContent(binx, biny); 
      }
      b[iy+1] = bicubicFunction( (x - h->GetXaxis()->GetBinCenter(binx0)) / xwidth,a); 
    }
    
    return bicubicFunction((y - h->GetYaxis()->GetBinCenter(biny0)) / ywidth,b);  


  }


  return 0; 

}

TH1 * UCorrelator::image::makePctileHist(const TH1 * h, const char * name , bool invert, int npctilebins, bool integral) 
{
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  fprintf(stderr,"makePctileHist requires at least ROOT 6.04\n"); 
  return 0; 
#else


  std::vector<double> vals; 
  vals.reserve(h->GetNbinsX()*h->GetNbinsY() * h->GetNbinsZ()); 

  double min = DBL_MAX; 
  double max = -DBL_MAX; 
  double sum = 0; 
  for (int i = 1; i <= h->GetNbinsX(); i++)
  {
    double iwidth = h->GetXaxis()->GetBinWidth(i);
    for (int j = 1; j <= h->GetNbinsY(); j++)
    {
      double jwidth = h->GetYaxis()->GetBinWidth(j);
      for (int k = 1; k <= h->GetNbinsZ(); k++)
      {
        double w = 1; 
        if (integral) 
        {
          double kwidth = h->GetZaxis()->GetBinWidth(k);
          w = iwidth*jwidth*kwidth; 
        }
        double val = w*h->GetBinContent(i,j,k); 
        if (val < min) min = val; 
        if (val > max) max = val; 
        sum += val; 
        vals.push_back(val); 
      }
    }
  }

  TH1D aux("aux_pctile_hist", "blah", npctilebins,min,max); 
  for (unsigned i = 0; i < vals.size(); i++) 
  {
    aux.Fill(vals[i], vals[i]/sum); 
  }

  TH1 * cumu = aux.GetCumulative(invert); 

  TString hist_name;
  if (name[0]=='_') hist_name.Form("%s%s", h->GetName(),name); 
  else hist_name = name;
  TString hist_title;
  hist_title.Form("Pctiles of %s", h->GetTitle());

  int ndims = h->ClassName()[2]-'0'; 

  TH1* hout = ndims == 1 ? (TH1*) new TH1F(hist_name.Data(), hist_title.Data(), h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()) : 
              ndims == 2 ? (TH1*) new TH2F(hist_name.Data(), hist_title.Data(), h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 
                                                                         h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax()) :
              ndims == 3 ? (TH1*) new TH3F(hist_name.Data(), hist_title.Data(), h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 
                                                                         h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax(),
                                                                         h->GetNbinsZ(), h->GetZaxis()->GetXmin(), h->GetZaxis()->GetXmax())
              : 0;

  //copy the actual axes to get the binning right!!! 
  //This also gets the axes names right
  h->GetXaxis()->Copy(*hout->GetXaxis()); 
  h->GetYaxis()->Copy(*hout->GetYaxis()); 
  h->GetZaxis()->Copy(*hout->GetZaxis()); 


  //copy the number of entries 
  hout->SetEntries(h->GetEntries()); 

  for (int i = 1; i <= h->GetNbinsX(); i++)
  {
    double iwidth = h->GetXaxis()->GetBinWidth(i);
    for (int j = 1; j <= h->GetNbinsY(); j++)
    {
      double jwidth = h->GetYaxis()->GetBinWidth(j);
      for (int k = 1; k <= h->GetNbinsZ(); k++)
      {
        double w = 1; 
        if (integral) 
        {
          double kwidth = h->GetZaxis()->GetBinWidth(k);
          w = iwidth*jwidth*kwidth; 
        }
        double val = w*h->GetBinContent(i,j,k); 
        double pctile = cumu->Interpolate(val); 
        if (invert) pctile = 1-pctile;
        hout->SetBinContent(i,j,k, pctile); 
      }
    }
  }

  delete cumu; 

  hout->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  if (ndims > 1) 
  {
    hout->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  }
  if (ndims > 2) 
  {
    hout->GetZaxis()->SetTitle("Percentile");
  }
  hout->SetMaximum(1); 
  hout->SetMinimum(0); 

  return hout; 
#endif

}

