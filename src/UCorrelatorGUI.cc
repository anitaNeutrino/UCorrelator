#include "UCorrelatorGUI.h" 
#include "TCanvas.h" 
#include "AnitaEventSummary.h" 
#include "WaveformCombiner.h" 
#include "TVirtualHistPainter.h"
#include "FilteredAnitaEvent.h" 


UCorrelator::gui::Map::Map(const TH2D & hist, const FilteredAnitaEvent * ev, WaveformCombiner * c, WaveformCombiner *cf, AnitaPol::AnitaPol_t pol) 
  : TH2D(hist), wfpad(0), f(ev), c(c), cf(cf), clicked(0), use_filtered(false),pol(pol)  
{
  SetStats(0);  
  SetDirectory(0); 

  /** figure out sun position */ 

  UsefulAdu5Pat pat(*ev->getGPS()); 

  double phi_sun,theta_sun; 
  pat.getSunPosition(phi_sun,theta_sun); 
  phi_sun *= 180/M_PI; 
  theta_sun *= 180/M_PI; 
  TMarker sun0(phi_sun,theta_sun,4); 
  TMarker sun1(phi_sun,theta_sun,7); 
  specials.push_back(sun0); 
  specials.push_back(sun1); 

  /** figure out if WAIS is in within ~ 1000 km */ 

  if (pat.getWaisDivideTriggerTimeNs() < 3.3e6)
  {
    double phi_wais, theta_wais; 
    pat.getThetaAndPhiWaveWaisDivide(phi_wais,theta_wais); 
    phi_wais *= 180/M_PI; 
    phi_wais *= 180/M_PI; 
    TMarker wais(phi_wais, theta_wais,29); 
    wais.SetMarkerColor(2); 
    specials.push_back(wais); 
  }

  /** figure out if LDB is in within ~ 1000 km */ 
  if (pat.getLDBTriggerTimeNs() < 3.3e6)
  {
    double phi_ldb, theta_ldb; 
    pat.getThetaAndPhiWaveLDB(phi_ldb,theta_ldb); 
    phi_ldb *= 180/M_PI; 
    phi_ldb *= 180/M_PI; 
    TMarker ldb(phi_ldb, theta_ldb,29); 
    ldb.SetMarkerColor(2); 
    specials.push_back(ldb); 
  }
}


void UCorrelator::gui::Map::addRough(const std::vector<std::pair<double,double> > & rough)
{
  for (unsigned i = 0; i < rough.size(); i++)
  {
    addRough( rough[i].first, rough[i].second); 
  }
}

void UCorrelator::gui::Map::addRough(double x, double y) 
{
  TMarker m(x, y,3); 
  m.SetMarkerStyle(2); 
  rough_m.push_back(m); 
}

void UCorrelator::gui::Map::addFine(const AnitaEventSummary::PointingHypothesis & p) 
{
  TMarker m(p.phi, -p.theta,2); 
  m.SetMarkerStyle(2); 
  m.SetMarkerColor(3); 
  fine_m.push_back(m); 

  double angle = 90. / TMath::Pi()* atan2(2*p.rho * p.sigma_theta * p.sigma_phi, p.sigma_phi * p.sigma_phi - p.sigma_theta * p.sigma_theta);
  TEllipse ell(p.phi, -p.theta, p.sigma_phi, p.sigma_theta, 0, 360, angle); 
  ell.SetFillStyle(0); 
  ell.SetLineColor(3); 
  fine_e.push_back(ell); 
}

UCorrelator::gui::Map::~Map()
{
  if (wfpad) delete wfpad; 
}


void UCorrelator::gui::Map::Paint(Option_t * opt) 
{
  //we ignore most of the options 
  GetPainter("colz2");
  fPainter->Paint("colz2"); 

  /* turn off peaks with np*/ 
  if (!strcasestr(opt,"np"))
  {
    for (size_t i = 0; i < rough_m.size(); i++) 
    {
      rough_m[i].SetMarkerSize(1+rough_m.size()-i); 
      rough_m[i].Draw(); 
    }

    for (size_t i = 0; i < fine_m.size(); i++) 
    {
      fine_m[i].Draw(); 
    }
  }

  /* turn off peaks with ne*/ 
  if (!strcasestr(opt,"ne"))
  {
    for (size_t i = 0; i < fine_e.size(); i++) 
    {
      fine_e[i].Draw(); 
    }
  }

  /** turn off specials with ns */ 

  if (!strcasestr(opt,"ns"))
  {
    for (size_t i = 0; i < specials.size(); i++) 
    {
      specials[i].Draw(); 
    }
  }
}

void UCorrelator::gui::Map::SetUseUnfiltered()
{
  use_filtered = false; 
}

void UCorrelator::gui::Map::SetUseFiltered()
{
  use_filtered = true; 

}

void UCorrelator::gui::Map::clear()
{
  rough_m.clear(); 
  fine_m.clear(); 
  fine_e.clear(); 
  specials.clear(); 
}

void UCorrelator::gui::Map::ExecuteEvent(int event, int px, int  py)
{
  if (event != kButton1Down && fPainter)
  {
    fPainter->ExecuteEvent(event,px,py); 
    return;
  }

  Double_t x  = gPad->PadtoX(gPad->AbsPixeltoX(px));
  Double_t y  = gPad->PadtoY(gPad->AbsPixeltoY(py));

  if (!wfpad) 
  {
    wfpad = new TCanvas(TString::Format("%s_zoom",GetName()),TString::Format("%s (clicked)", GetTitle()), 800,800); 
  }

  wfpad->Clear(); 
  WaveformCombiner * wf = use_filtered ? cf : c; 
  wf->combine(x, -y, f, (AnitaPol::AnitaPol_t) pol, 0); 


  bool do_deconvolve = wf->getDeconvolved(); 

  wfpad->Divide(2,do_deconvolve?2:1); 

  wfpad->cd(1); 
  wf->getCoherent()->drawEven(); 
  wfpad->cd(2); 
  wf->getCoherent()->drawPowerdB(); 

  if (do_deconvolve)
  {
    wfpad->cd(3); 
    wf->getDeconvolved()->drawEven(); 
    wfpad->cd(4); 
    wf->getDeconvolved()->drawPowerdB(); 
  }
}
