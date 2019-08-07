#include "EASFitter.h" 
#include "FFTtools.h" 
#include "DeltaT.h" 
#include "AntennaPositions.h" 
#include "TH2.h" 
#include "AnalysisWaveform.h" 
#include "TGraph2D.h" 
#include "TCanvas.h" 
#include "TFile.h" 
#include "Minuit2/Minuit2Minimizer.h" 
#include "FilteredAnitaEvent.h" 
#include "FreqDomainFunction.h" 
#include "TimeDependentAverage.h" 
#include "RawAnitaHeader.h" 
#include "ResponseManager.h" 



static TH2D * hpower =0;
static TH2D * hphase =0;
static std::vector<TGraph*> vpower; 
static std::vector<TGraph*> vphase; 


static double getOffAxisPowerMultiplier(int pol, int plane, double f, double deg)
{

  static TFile foff(Form("%s/share/data/A3offaxis.root", getenv("ANITA_UTIL_INSTALL_DIR"))); 
  static bool init = false; 
  static TGraph2D * data[2][2]; 

  if (!init) 
  {
    data[0][1]= (TGraph2D*) foff.Get("Hpol_Hplane");
    data[0][0]= (TGraph2D*) foff.Get("Hpol_Eplane");
    data[1][1]= (TGraph2D*) foff.Get("Vpol_Hplane");
    data[1][0] = (TGraph2D*) foff.Get("Vpol_Eplane");
    init = true; 
  }

  return pow(10,data[pol][plane]->Interpolate(f, deg)/10);
}

static void setupHists()
{

  TFile f(Form("%s/share/data/templates/crTmpltsA3.root",getenv("ANITA_UTIL_INSTALL_DIR"))); 

  for (int i = 0; i < 27; i++) 
  {
    TGraph * g = FFTtools::cropWave((TGraph*) f.Get(Form("efield%d",i)),0,600); 
    g->SetTitle(Form("Template %d",i)); 

    AnalysisWaveform * wf = new AnalysisWaveform(g->GetN(),g->GetY(), g->GetX()[1]-g->GetX()[0], g->GetX()[0]); 
    wf->setTitle(g->GetTitle()); 

    int roll; 
    wf->hilbertEnvelope()->peakVal(&roll);  
    FFTtools::rotate(wf->updateEven(), -roll); 

    TGraphAligned * gpower = (TGraphAligned*) wf->power(); 
    TGraphAligned * gphase = new TGraphAligned(*wf->phase()); 

    //find the peak of the power spectrum 
    int peak; 
    double val = gpower->peakVal(&peak); 
    FFTtools::unwrap(gphase->GetN(), gphase->GetY(), 2*TMath::Pi()); 
    for (int j = 0; j < gpower->GetN(); j++) 
    {
      if (j < peak) continue; 

      if (gpower->GetY()[j] / val < 1e-3) 
      {
        gphase->Set(j); 
        break; 
      }
    }

    gpower->SetTitle(g->GetTitle());
    gphase->SetTitle(g->GetTitle());
    vpower.push_back(gpower);
    vphase.push_back(gphase); 
  }

  double df = vpower[0]->GetX()[1] - vpower[0]->GetX()[0];

  hpower = new TH2D("hpower","POWER;freq;T", vpower[0]->GetN(), -df/2, vpower[0]->GetX()[vpower[0]->GetN()-1]+df/2,
      27,-13.5,13.5); 
  hphase = new TH2D("hphase","PHASE;freq;T", vpower[0]->GetN(), -df/2, vpower[0]->GetX()[vpower[0]->GetN()-1]+df/2,
      27,-13.5,13.5); 

  hpower->SetDirectory(0); 
  hphase->SetDirectory(0); 

  for (int x = 1; x <= hpower->GetNbinsX(); x++) 
  {
    for (int y = 1; y <= hpower->GetNbinsY(); y++)
    {
       hpower->SetBinContent(x,y, vpower[y-1]->GetY()[x-1]); 
       hphase->SetBinContent(x,y, vphase[y-1]->GetN() < x ? vphase[y-1]->GetY()[vphase[y-1]->GetN()-1] : vphase[y-1]->GetY()[x-1]); 
    }
  }

}

static std::complex<double> crfitfn (double f, const double * pars) 
{

  if (!hpower) setupHists(); 

  double T = *pars; 
  if (T < -13 || T > 13) return 0; 
  if (f < 0 || f> hpower->GetXaxis()->GetBinCenter(hpower->GetNbinsX())) return 0; 
  double dphi = pars[1];
  double dtheta = pars[2];
  int ipol = int(pars[3]); 
  //printf("%g %g %g %d\n", T,dphi,dtheta,ipol); 

  double atten_phi = dphi == 0 ? 1 :  getOffAxisPowerMultiplier(ipol, ipol == 0,f, dphi); 
  double atten_theta = dtheta == 0 ? 1 : getOffAxisPowerMultiplier(ipol, ipol == 1,f, dtheta); 

//  printf(" %g %g\n",f,T); 
  double pow = hpower->Interpolate(f,T) * atten_phi*atten_theta; 
  double phase = hphase->Interpolate(f,T); 
  double amp = sqrt(pow); 

  return std::complex<double>( cos(phase) * amp, sin(phase)*amp); 
}




class EASFitFn : public ROOT::Math::IBaseFunctionMultiDim
{


  public:
  virtual ~EASFitFn() 
  {

    for (unsigned i = 0; i < fns.size(); i++)
    {
      delete fns[i]; 
      delete wfs[i]; 
    }
  }

  
  EASFitFn(int nants, const int * ants, const FilteredAnitaEvent * ev, int pol, double t, double phi, double theta, const AnitaResponse::ResponseManager  * rm, bool freqonly = false, bool dedisperse = false) 
    : dedisperse(dedisperse), freqonly(freqonly), nants(nants), pol(pol), ants(ants), rms(nants), event(ev),wfs(nants),  fns(nants), t(t), phi(phi), theta(theta), rm(rm)
  {

    AnitaResponse::AllPassDeconvolution allpass;
    for (int i = 0; i < nants; i++) 
    {
      rms[i] = freqonly ? 0 :  UCorrelator::TimeDependentAverageLoader::getRMS(t,0,ants[i]); 
      wfs[i] = new AnalysisWaveform(*ev->getRawGraph(ants[i],AnitaPol::kHorizontal));
      if (!freqonly && dedisperse) rm->response(0,ants[i])->deconvolveInPlace(wfs[i], &allpass); 
      fns[i] = new FreqDomainFunction(crfitfn, 4, rm->response(0,ants[i]),wfs[i]->deltaT()/4,-wfs[i]->even()->GetX()[wfs[i]->Neven()-1], wfs[i]->even()->GetX()[wfs[i]->Neven()-1]);
      if (!freqonly && dedisperse) fns[i]->dedisperseResponse(true); 
    }
  }
  virtual ROOT::Math::IBaseFunctionMultiDim * Clone() const { return new EASFitFn(rms.size(), ants, event, pol, t, phi,theta,rm,freqonly,dedisperse); } 
  virtual unsigned int NDim() const { return 2 *nants+3; } 
  

  virtual double DoEval(const double *x) const
  {
    double A = x[0]; 
    double t0 = x[1]; 
    double T = x[2]; 
    double chisq = 0; 
//    printf("T: %g:\n", T); 

    for (int i = 0; i < nants; i++) 
    {
      double delta_t = UCorrelator::getDeltaT(ants[0], ants[i], phi, theta, AnitaPol::kHorizontal,true); 
      double dphi = FFTtools::wrap(UCorrelator::AntennaPositions::instance()->phiAnt[0][ants[i]] -phi,360,0); 
      double dtheta = 10 - theta;

      double p[6] = { A*x[2*i+3], freqonly ? 0 : t0-delta_t + x[2*i+4],  T, dphi, dtheta,(double)pol}; 
//      printf("%g %g %g %g %g\n", p[0],p[1], p[2],p[3],p[4]); 
      fns[i]->setParameters(p); 

      if (freqonly) 
      {
        double df = wfs[i]->deltaF(); 
        for (int j = 0; j < wfs[i]->Nfreq(); j++) 
        {
            chisq += pow(wfs[i]->power()->GetY()[j] - fns[i]->evalPower(j*df),2); 
        }
      }
      else
      {
        for (int j = 0; j < wfs[i]->Neven(); j++) 
        {
          if (dedisperse) 
            chisq += pow(wfs[i]->even()->GetY()[j]  - fns[i]->eval(wfs[i]->even()->GetX()[j]),2)/((wfs[i]->Neven() - 3) *rms[i]*rms[i]);
          else 
            chisq += pow(wfs[i]->uneven()->GetY()[j]  - fns[i]->eval(wfs[i]->uneven()->GetX()[j]),2)/((wfs[i]->Nuneven() - 3) *rms[i]*rms[i]);
        }
      }
    }
    return chisq; 
  }


  void plot(TCanvas *c, std::vector<TObject*> *save = 0)
  {
    c->Divide(nants/3,3); 

    for (int i = 0; i < nants; i++) 
    {
      c->cd(i+1); 
      AnalysisWaveform *wf =   wfs[i]; 
      wf->updateEven()->SetTitle(Form("Ant %d", ants[i])); 
      wf->updateEven()->setPlottingLimits(1.1,true,20); 
      if (save) 
      {
        TGraphAligned * g = new TGraphAligned(*wf->even());
        save->push_back(g); 
        g->Draw("al"); 
      }
      else
      {
        wf->drawEven("al"); 
      }
      TF1 * f = fns[i]->makeTF1(Form("fn_%d",i)); 
      if (save) save->push_back(f); 
      f->Draw("same"); 
    }
  }



  bool dedisperse; 
  bool freqonly; 
  int nants; 
  int pol; 
  const int* ants; 
  std::vector<double> rms; 
  const FilteredAnitaEvent * event; 
  std::vector<AnalysisWaveform *> wfs; 
  mutable std::vector<FreqDomainFunction *> fns; 
  double t;
  double phi,theta;
  const AnitaResponse::ResponseManager * rm; 

};



int UCorrelator::EASFitter::fitEvent(int nants, const int * ants, 
                                    const FilteredAnitaEvent * event, 
                                    double phi, double theta, bool dedisperse) 
{

  EASFitFn freq_fit(nants, ants, event, AnitaPol::kHorizontal, event->getHeader()->triggerTime, phi,theta,rm,true,dedisperse); 
  EASFitFn time_fit(nants, ants, event, AnitaPol::kHorizontal, event->getHeader()->triggerTime, phi,theta,rm,false,dedisperse); 

  double freq_scale = freq_fit.fns[0]->waveform()->power()->peakVal(); 

  int toffset_loc;
  double scale = time_fit.fns[0]->waveform()->even()->peakVal(&toffset_loc,0,-1,true); 

  for (int side_of_cone = -1; side_of_cone <=1; side_of_cone+=2) 
  {

    //start with a frequency domain fit, which doesn't care about polarity 
    ROOT::Minuit2::Minuit2Minimizer min_freq; 
    AnalysisWaveform * wf = (AnalysisWaveform*) freq_fit.wfs[0]->freq(); 
    int loc; 
    double fA = wf->even()->peakVal(&loc, 0,-1,true) / freq_scale ; 
    int ivar = 0;
    min_freq.SetLimitedVariable(ivar++, "A", fA, fA*0.1, 0, 2*fA);
    min_freq.SetFixedVariable(ivar++, "t0", 0); //doesn't matter for frequency
    min_freq.SetLimitedVariable(ivar++,"T", side_of_cone*7,0.1,side_of_cone > 0 ? 0 : -13, side_of_cone > 0 ? 13 : 0); 
      
    for (int i =0; i < nants; i++) 
    {
    //    min.SetLimitedVariable(ivar++, Form("relA_%d",i), 1,0.1, 0.8, 1.2); 
        min_freq.SetFixedVariable(ivar++, Form("relA_%d",i), 1); 
        min_freq.SetFixedVariable(ivar++, Form("tadj_%d",i), 0); //doesn't matter for frequency fit
    }

    min_freq.SetFunction(freq_fit); 
    min_freq.Minimize(); 

    printf("Side of Cone: %d, status: %d.  T: %g, A: %g\n",
        side_of_cone, min_freq.Status(), min_freq.X()[2], min_freq.X()[0]); 


    for (int polarity = -1; polarity <=1; polarity+=2) 
    {


      ROOT::Minuit2::Minuit2Minimizer min_time; 
      
      min_time.SetFunction(time_fit); 


      const AnalysisWaveform * wf = (AnalysisWaveform*) time_fit.wfs[0]; 
      int loc; 
      double A = wf->even()->peakVal(&loc, 0,-1,true) / scale ; 
      double t0 = wf->even()->GetX()[loc]; 

      int ivar = 0;
      min_time.SetLimitedVariable(ivar++, "A", A, A*0.1, 0, 2*A);
      min_time.SetVariable(ivar++, "t0", t0, 0.5); 
      min_time.SetLimitedVariable(ivar++,"T", 7,0.1,-13,13); 
      
      for (int i =0; i < nants; i++) 
      {
        printf("ant: %d\n", ants[i]); 
    //    min.SetLimitedVariable(ivar++, Form("relA_%d",i), 1,0.1, 0.8, 1.2); 
        min_time.SetFixedVariable(ivar++, Form("relA_%d",i), 1); 
        min_time.SetLimitedVariable(ivar++, Form("tadj_%d",i), 0, 0.1 , -1, 1); 
      }

      min_time.Minimize(); 
      min_time.PrintResults(); 

    }
  }
  return 0; 
}


UCorrelator::EASFitter::~EASFitter() 
{
  for (unsigned i = 0; i < plots.size(); i++) delete plots[i]; 
  for (unsigned i = 0; i < save.size(); i++) delete plots[i]; 
}



  
