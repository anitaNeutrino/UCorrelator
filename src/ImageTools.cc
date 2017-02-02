#include "TH2.h" 
#include "UCImageTools.h" 
#include "TString.h" 
#include <algorithm>



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
    case 'L': 
      return new TH1L(name,title,nbins,xmin,xmax); 
    case 'F': 
      return new TH1F(name,title,nbins,xmin,xmax); 
    case 'D': 
    default:
      return new TH1D(name,title,nbins,xmin,xmax); 
  }
}



TH1* UCorrelator::image::getPctileProjection(const TH2 * H, int axis, double pct) 
{
  TString name; 
  TString title; 
  name.Form("%s_%g_proj_%d", H->GetName(), pct, axis); 
  title.Form("%s %g%% projection of axis %d", H->GetTitle(), 100*pct, axis); 
  
  TH1 * h = makeTH1(H->GetClassName()[3], name.Data(), title.Data(); H->GetAxis(axis)->GetNbins(), H->GetAxis(axis)->GetXmin(), H->GetAxis(axis)->GetXmax()); 

  int nOrth = h->GetXaxis(3-axis)->GetNbins(); 

  for (int i = 1; i <= h->GetNbinsX(); i++) 
  {
    std::vector<double> v(nOrth); 
    for (int j = 1; j <= nOrth; j++)
    {
      v[j-1] = H->GetBinContent( axis == 1 ? i : j, axis == 1 ? j : i); 
    }
    std::nth_element(v.begin(), v.begin() + (nOrth-1) * pct, v.end()); 
    h->SetBinContent(i, v[(nOrth-1)*pct]; 
  }

  return h; 
} 

