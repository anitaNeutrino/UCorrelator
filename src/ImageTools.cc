#include "TProfile2D.h"
#include "TH2.h" 
#include "UCImageTools.h" 
#include "TString.h" 
#include <algorithm>
#include <set>
#include <math.h>



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




double UCorrelator::getZRMS(const TH2* hist) {

  int numBinsX = hist->GetNbinsX();
  int numBinsY = hist->GetNbinsY();

  double sum2 = 0;

  for (int binX=0; binX<numBinsX; binX++) {
    for (int binY=0; binY<numBinsY; binY++) {
      double binValue = double(hist->GetBinContent(binX,binY));
      sum2 += pow(binValue,2);
    }
  }

  return sqrt(sum2/(numBinsX*numBinsY));

}


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

