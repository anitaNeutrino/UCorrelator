#include "TH2.h" 
#include "UCImageTools.h" 
#include "TString.h" 
#include <algorithm>
#include <set>



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


TH1* UCorrelator::image::getPctileProjection(const TH2 * H, int axis, double pct, bool ignoreEmpty) 
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

