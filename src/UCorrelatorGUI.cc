#include "UCorrelatorGUI.h" 
#include "TCanvas.h" 
#include "AnitaEventSummary.h" 
#include "FFTtools.h" 
#include "AnalysisWaveform.h" 
#include "Analyzer.h" 
#include "WaveformCombiner.h" 
#include "TVirtualHistPainter.h"
#include "FilteredAnitaEvent.h" 


ClassImp(UCorrelator::gui::Map); 
ClassImp(UCorrelator::gui::SummaryText); 



UCorrelator::gui::Map::Map(const TH2D & hist, const FilteredAnitaEvent * ev, WaveformCombiner * c, WaveformCombiner *cf, AnitaPol::AnitaPol_t pol, const AnitaEventSummary* sum ) 
  : TH2D(hist), wfpad(0), f(ev), c(c), cf(cf), clicked(0), use_filtered(false),pol(pol)  , heading_axis(GetXaxis()->GetXmin(), GetYaxis()->GetXmax(), GetXaxis()->GetXmax(), GetYaxis()->GetXmax(), GetXaxis()->GetXmin()-ev->getGPS()->heading, GetXaxis()->GetXmax() - ev->getGPS()->heading,510,"=")
{
  SetStats(0);  
  SetDirectory(0); 

  /** figure out sun position */ 



  double phi_sun = sum->sun.phi;
  double theta_sun = -sum->sun.theta; 
  phi_sun = FFTtools::wrap(phi_sun,360); 
  specials.push_back(TMarker(phi_sun,theta_sun,4)); 
  specials.push_back(TMarker(phi_sun,theta_sun,7)); 


  if (sum->mc.phi >-999)
  {
//    printf("MC position: phi=%g theta=%g\n",sum->mc.phi, sum->mc.theta ); 
    TMarker mc(sum->mc.phi, -sum->mc.theta,29); 
    mc.SetMarkerColor(2); 
    specials.push_back(mc); 
  }
  else
  {
    /** figure out if WAIS is (roughly) above horizon*/ 
    if (sum->wais.theta < 8)
    {
      TMarker wais(sum->wais.phi, -sum->wais.theta,29); 
      wais.SetMarkerColor(2); 
      specials.push_back(wais); 
    }

    /** figure out if LDB is (roughly) above horizon */
    if (sum->ldb.theta < 8)
    {
      TMarker ldb(sum->ldb.phi, -sum->ldb.theta,29); 
      ldb.SetMarkerColor(2); 
      specials.push_back(ldb); 
    }
  }
  heading_axis.SetTitle("heading"); 

  coherent = 0; 
  deconvolved = 0; 
}

UCorrelator::gui::Map::Map(const Map & other) 
  : TH2D(other) , rough_m(other.rough_m), fine_m(other.fine_m), fine_e(other.fine_e), 
   specials(other.specials),  wfpad(other.wfpad), f(other.f), c(other.c), cf(other.cf),
   clicked(other.clicked), use_filtered(other.use_filtered),  pol(other.pol), last_theta(other.last_theta),
   last_phi(other.last_phi), heading_axis(GetXaxis()->GetXmin(), GetYaxis()->GetXmax(), GetXaxis()->GetXmax(), GetYaxis()->GetXmax(), GetXaxis()->GetXmin()-f->getGPS()->heading, GetXaxis()->GetXmax() - f->getGPS()->heading,510,"=")
{
  SetStats(0);  
  SetDirectory(0); 
  heading_axis.SetTitle("heading"); 
  coherent = other.coherent ? new AnalysisWaveform(*other.coherent) : 0; 
  deconvolved = other.deconvolved ? new AnalysisWaveform(*other.deconvolved) : 0; 

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
  if (clicked) delete clicked; 
  if (coherent) delete coherent; 
  if (deconvolved) delete deconvolved; 
}


void UCorrelator::gui::Map::Paint(Option_t * opt) 
{
  //we ignore most of the options 
  GetPainter(opt);
  fPainter->Paint(opt); 
  if (!strcasestr(opt,"nh"))
  {
    heading_axis.Draw(); 
  }

  /* turn off peaks with np*/ 
  if (!strcasestr(opt,"np"))
  {
    for (size_t i = 0; i < rough_m.size(); i++) 
    {
      if (rough_m[i].GetX() < GetXaxis()->GetXmin()) continue;
      if (rough_m[i].GetX() > GetXaxis()->GetXmax()) continue;
      if (rough_m[i].GetY() < GetYaxis()->GetXmin()) continue;
      if (rough_m[i].GetY() > GetXaxis()->GetXmax()) continue;

      rough_m[i].SetMarkerSize(1+rough_m.size()-i); 
      rough_m[i].Draw(); 
    }

    for (size_t i = 0; i < fine_m.size(); i++) 
    {
      if (fine_m[i].GetX() < GetXaxis()->GetXmin()) continue;
      if (fine_m[i].GetX() > GetXaxis()->GetXmax()) continue;
      if (fine_m[i].GetY() < GetYaxis()->GetXmin()) continue;
      if (fine_m[i].GetY() > GetXaxis()->GetXmax()) continue;

 
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
      if (specials[i].GetX() < GetXaxis()->GetXmin()) continue;
      if (specials[i].GetX() > GetXaxis()->GetXmax()) continue;
      if (specials[i].GetY() < GetYaxis()->GetXmin()) continue;
      if (specials[i].GetY() > GetXaxis()->GetXmax()) continue;

 
      specials[i].Draw(); 
    }
  }

  if (clicked) clicked->Draw(); 
}

void UCorrelator::gui::Map::SetUseUnfiltered()
{
  use_filtered = false; 
  if (wfpad) drawWf(last_phi, last_theta); 
}

void UCorrelator::gui::Map::SetUseFiltered()
{
  use_filtered = true; 
  if (wfpad) drawWf(last_phi, last_theta); 

}

void UCorrelator::gui::Map::clear()
{
  rough_m.clear(); 
  fine_m.clear(); 
  fine_e.clear(); 
  specials.clear(); 
}

void UCorrelator::gui::Map::closeCanvas()
{
  TMarker * deleteme = clicked; 
  clicked = 0;
  wfpad = 0; 
  delete deleteme; 
}

void UCorrelator::gui::Map::ExecuteEvent(int event, int px, int  py)
{
  if (event != kButton1Down && fPainter)
  {
//    fPainter->ExecuteEvent(event,px,py); //??? 
    return;
  }

  Double_t x  = gPad->PadtoX(gPad->AbsPixeltoX(px));
  Double_t y  = gPad->PadtoY(gPad->AbsPixeltoY(py));

  TVirtualPad * ours = gPad; 
  drawWf(x,-y); 

  if (clicked) delete clicked; 
  clicked = new TMarker(x,y,5); 
  ours->Modified(); 
  ours->Update(); 

}

void UCorrelator::gui::Map::drawWf(double phi, double theta) 
{
  last_theta=theta;
  last_phi = phi; 
  if (!wfpad || !wfpad->GetCanvasImp()) 
  {
    if (!wfpad) delete wfpad; 
    wfpad = new TCanvas(TString::Format("%s_zoom",GetName()),TString::Format("%s (clicked)", GetTitle()), 800,800); 
    wfpad->Connect("Closed()", "UCorrelator::gui::Map",this,"closeCanvas()"); 
  }

  wfpad->Clear(); 
  WaveformCombiner * wf = use_filtered ? cf : c; 
  wf->combine(phi, theta, f, (AnitaPol::AnitaPol_t) pol, 0); 


  if (coherent) delete coherent; 
  if (deconvolved) delete deconvolved; 
  coherent = new AnalysisWaveform(*wf->getCoherent()); 
  deconvolved = wf->getDeconvolved() ? new AnalysisWaveform(*wf->getDeconvolved()) : 0; 

  coherent->updateEven()->SetTitle(TString::Format("Coherent (#phi=%g, #theta=%g) (%s)",phi,theta, use_filtered ? "filtered" :"unfiltered")); 

  wfpad->Divide(2,deconvolved?2:1); 
  wfpad->cd(1); 
  coherent->drawEven(); 
  wfpad->cd(2); 
  coherent->drawPowerdB(); 

  if (deconvolved)
  {
    deconvolved->updateEven()->SetTitle(TString::Format("Deconvolved (#phi=%g, #theta=%g) (%s)",phi,theta, use_filtered ? "filtered" :"unfiltered")); 
    wfpad->cd(3); 
    deconvolved->drawEven(); 
    wfpad->cd(4); 
    deconvolved->drawPowerdB(); 
  }


  wfpad->Update(); 
}


UCorrelator::gui::SummaryText::SummaryText(int i,AnitaPol::AnitaPol_t pol, const Analyzer *analyzer, int use_filtered)
  : TPaveText(i/double(analyzer->getSummary()->nPeaks[(int)pol]), 0, (i+1) /double(analyzer->getSummary()->nPeaks[(int) pol]),1)

{
  int ipol = (int) pol; 
  const AnitaEventSummary * ev = analyzer->getSummary(); 
  double rough_theta = analyzer->getRoughTheta(pol,i); 
  double rough_phi = analyzer->getRoughPhi(pol,i); 
  AddText(TString::Format("#phi: %0.2f (rough) , %0.3f (fine)", rough_phi, ev->peak[ipol][i].phi)); 
  AddText(TString::Format("#theta: %0.2f (rough) , %0.3f (fine)", rough_theta, ev->peak[ipol][i].theta)); 
  AddText(TString::Format("peak val: %f", ev->peak[ipol][i].value)); 
  if(use_filtered == 0)
  {
    AddText(TString::Format("peak_{hilbert}:  %0.3f (coher), %0.3f (deconv)", ev->coherent[ipol][i].peakHilbert, ev->deconvolved[ipol][i].peakHilbert)); 
    AddText(TString::Format("stokes: (coher): (%0.3g, %0.3g, %0.3g, %0.3g)", ev->coherent[ipol][i].I, ev->coherent[ipol][i].Q, ev->coherent[ipol][i].U, ev->coherent[ipol][i].V));
    AddText(TString::Format("stokes: (deconv): (%0.3g, %0.3g, %0.3g, %0.3g)", ev->deconvolved[ipol][i].I, ev->deconvolved[ipol][i].Q, ev->deconvolved[ipol][i].U, ev->deconvolved[ipol][i].V));
    AddText(TString::Format("impulsivity measure: %0.3g (coher), %0.3g (deconv)", ev->coherent[ipol][i].impulsivityMeasure, ev->deconvolved[ipol][i].impulsivityMeasure));
    AddText(TString::Format("SNR %0.3g (coher), %0.3g (deconv)", ev->coherent[ipol][i].snr, ev->deconvolved[ipol][i].snr));
  }
  else 
  {
    AddText(TString::Format("peak_{hilbert}:  %0.3f (coher), %0.3f (deconv)", ev->coherent_filtered[ipol][i].peakHilbert, ev->deconvolved_filtered[ipol][i].peakHilbert)); 
    AddText(TString::Format("max stokes: (coher): (%0.3g, %0.3g, %0.3g, %0.3g)", ev->coherent_filtered[ipol][i].max_dI, ev->coherent_filtered[ipol][i].max_dQ, ev->coherent_filtered[ipol][i].max_dU, ev->coherent_filtered[ipol][i].max_dV));
    AddText(TString::Format("max stokes: (deconv): (%0.3g, %0.3g, %0.3g, %0.3g)", ev->deconvolved_filtered[ipol][i].max_dI, ev->deconvolved_filtered[ipol][i].max_dQ, ev->deconvolved_filtered[ipol][i].max_dU, ev->deconvolved_filtered[ipol][i].max_dV));
    AddText(TString::Format("impulsivity measure: %0.3g (coher), %0.3g (deconv)", ev->coherent_filtered[ipol][i].impulsivityMeasure, ev->deconvolved_filtered[ipol][i].impulsivityMeasure));
    AddText(TString::Format("SNR %0.3g (coher), %0.3g (deconv)", ev->coherent_filtered[ipol][i].snr, ev->deconvolved_filtered[ipol][i].snr));
  }
  AddText(TString::Format("position: %0.3f N, %0.3f E, %0.3f m", ev->peak[ipol][i].latitude, ev->peak[ipol][i].longitude, ev->peak[ipol][i].altitude)); 
}
