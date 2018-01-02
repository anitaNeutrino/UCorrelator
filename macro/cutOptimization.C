

#include "macro/eff2d.C" 



static AnitaTemplateMachine * m = 0;
static UCorrelator::Analyzer * a = 0;

UCorrelator::AnalysisConfig cfg; 
FilterStrategy * strategy = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2"); 


std::map<int,double> crays; 

double getCray(int run, int event, int pk = 0)
{
  if (!m) 
  {
    m = new AnitaTemplateMachine; 
    m->loadTemplates(); 
    cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
    cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
    cfg.combine_nantennas = 15; 
    cfg.zoomed_nant = 15; 
    a = new UCorrelator::Analyzer(&cfg,true); 
  }

  if (crays.count(event)) return crays[event]; 

  AnitaDataset d(run);
  d.getEvent(event); 

  FilteredAnitaEvent ev(d.useful(), strategy, d.gps(), d.header()); 

  AnitaEventSummary sum; 
  a->analyze(&ev,&sum); 

  AnitaTemplateSummary ats; 
  m->doTemplateAnalysis(a->getCoherent(AnitaPol::kHorizontal,pk),0,pk,&ats); 

  double ans =  ats.coherent[0][pk].cRay[4]; 
  crays[event] = ans; 
  return ans; 
}



void makeHists(const char * output = "cutopt.root", bool make_thermal_hists = false, bool make_cray = false)
{

  TFile of(output,"RECREATE"); 

  //Get the efficiency as a function of O/F  
  TChain cmc("overlap"); 
  cmc.Add("source_maps_eval_nosubtract/*sim.root"); 


  //we need the total weight, we can get that from the thermal tree
  //just use the first 100 MC runs for now
  TChain cmc_th("simulation"); 
  cmc_th.Add("thermalTrees/simulated_1-100*.root"); 

  TH1 * wgtHist = new TH1D("mcwgthist","mc weight hist", 100,0,1); 
  cmc_th.Draw("weight >> mcwgthist","weight * isMostImpulsive","goff"); 
  printf("runs 1-100 total weight: %g\n",wgtHist->Integral()); 
  wgtHist->Write(); 


  TCanvas * cefficiency = new TCanvas("cefficiency","Efficiency",800,600); 
  TH2 * efficiency = eff2D(&cmc,"-log(O+exp(-40))", 40,0, 40, "F",50,2,10, wgtHist->Integral(), "weight * (O>=0) * (run <=100)"); 
  efficiency->SetStats(false); 
  efficiency->SetName("eff"); 
  efficiency->Write(); 
  efficiency->DrawCopy("colz"); 





  TChain * cdata = new TChain("overlap"); 
  cdata->Add("source_maps_eval_nosubtract/*full.root"); 


  TCanvas * canthro = new TCanvas("canthro","canthro",800,600); 
  canthro->Divide(3,1); 

  TH2 * hvanthro = new TH2D("hvanthro","H+V", 70,-20,50,150,2,10); 
  hvanthro->GetYaxis()->SetTitle("F"); 
  hvanthro->GetXaxis()->SetTitle("-log(O + eps)"); 
  cdata->Draw("TMath::Min(F,9.99):TMath::Min(-log( exp(-50) + abs(O)),49.99) >> hvanthro", "(O >-1) * (O >1 || F <2.8) * (nclustered[2] < 100)","goff");  

  int n_pass_thermal_hv_highoverlap = cdata->Draw("O", "(O >-1) * (O >1)  && (F > 3.2) * (nclustered[2] < 100)","goff");  
  int n_pass2_thermal_hv_highoverlap = cdata->Draw("O", "(O >-1) * (O >1)  && (F > 2.8) * (nclustered[2] < 100)","goff");  
  int n_fail_thermal_hv_highoverlap = cdata->Draw("O", "(O >-1) * (O >1)  && (F < 2.8) * (nclustered[2] < 100)","goff");  
  printf("HV tau ~= %d/(%d  - %d) = %g-%g\n", n_fail_thermal_hv_highoverlap, n_pass2_thermal_hv_highoverlap,  n_pass_thermal_hv_highoverlap, double(n_fail_thermal_hv_highoverlap) / n_pass_thermal_hv_highoverlap, double(n_fail_thermal_hv_highoverlap) / n_pass2_thermal_hv_highoverlap); 


  TH2 * hanthro = new TH2D("hanthro","H", 70,-20,50,150,2,10); 
  hanthro->GetYaxis()->SetTitle("F"); 
  hanthro->GetXaxis()->SetTitle("-log(O + eps)"); 
  cdata->Draw("TMath::Min(F,9.99):TMath::Min(-log( exp(-50) + abs(O)),49.99) >> hanthro", "(O >-1) * (O >1 || F <2.8) * (pol == 0) * (nclustered[2] < 100)","goff");  


  TH2 * vanthro = new TH2D("vanthro","V", 70,-20,50,150,2,10); 
  vanthro->GetYaxis()->SetTitle("F"); 
  vanthro->GetXaxis()->SetTitle("-log(O + eps)"); 
  cdata->Draw("TMath::Min(F,9.99):TMath::Min(-log( exp(-50) + abs(O)),49.99) >> vanthro", "(O >-1) * (O >1 || F <2.8) * (pol == 1) * (nclustered[2] < 100)","goff");  

  int n_pass_thermal_v_highoverlap = cdata->Draw("O", "(O >-1) * (O >1)  && (F > 3.2) * (nclustered[2] < 100) && pol == 1","goff");  
  int n_fail_thermal_v_highoverlap = cdata->Draw("O", "(O >-1) * (O >1)  && (F < 2.8) * (nclustered[2] < 100) && pol","goff");  
  printf("V tau ~= %d/%d = %g\n", n_fail_thermal_v_highoverlap, n_pass_thermal_v_highoverlap, double(n_fail_thermal_v_highoverlap) / n_pass_thermal_v_highoverlap); 



  hvanthro->Write(); 
  hanthro->Write(); 
  vanthro->Write(); 
  hvanthro->SetStats(false); 
  hanthro->SetStats(false); 
  vanthro->SetStats(false); 



  TPaveText * hidden = new TPaveText(0, 2.8,50,10); 
  hidden->AddText("hidden"); 
  canthro->cd(1); 
  hvanthro->DrawCopy("colz"); 
  hidden->Draw("same"); 

  canthro->cd(2); 
  hanthro->DrawCopy("colz"); 
  hidden->Draw("same"); 

  canthro->cd(3); 
  vanthro->DrawCopy("colz"); 
  hidden->Draw("same"); 


  if (make_cray) 
  {

    cdata->SetScanField(0); 


    TCanvas * chpolall = new TCanvas("chpol_all","chpol_all",1024,800); 
    cdata->SetMarkerStyle(21); 
    int npoints = cdata->Draw("F: (O <=0) * (50) +  (O>0) * -log( exp(-50) +O):getCray(run,event)", "(O >-1) * (pol == 0)  * (nclustered[2] < 100)","colz");  

    TGraph2D * gdata = new TGraph2D(npoints); 

    for (int i = 0; i < npoints; i++)
    {
      gdata->SetPoint(i, cdata->GetV1()[i], cdata->GetV2()[i], cdata->GetV3()[i]); 
    }

    gdata->GetXaxis()->SetTitle("-log O"); 
    gdata->GetYaxis()->SetTitle("F"); 
    gdata->GetZaxis()->SetTitle("CR template 4 XCorr"); 
    gdata->SetName("all_hpol"); 
    gdata->SetTitle("HPol Events"); 
    gdata->Write("all_hpol"); 

    cdata->Scan("run:event:F:O:theta:getCray(run,event)", "(O >-1) * (O < 1) * (pol == 0) * (nclustered[2] < 100) ","goff");  
  }

  if (make_thermal_hists)
  {
    // for setting thermal cut 
    TChain cdata_th("anita3"); 
    cdata_th.SetProof(1); 
    cdata_th.Add("thermalTrees/a3all*.root"); 
    TH1 * thermal_hv = new TH1D("thermal_hv","Thermal Cut (above horizon sideband, pass quality cut, H+V)",500,-15,5); 
    thermal_hv->GetXaxis()->SetTitle("F"); 

    TH1 * thermal_h = new TH1D("thermal_h","Thermal Cut (above horizon sideband, pass quality cut H)",500,-15,5); 
    thermal_h->SetLineColor(11); 
    thermal_h->GetXaxis()->SetTitle("F"); 

    TH1 * thermal_v = new TH1D("thermal_v","Thermal Cut (above horizon sideband, pass quality cut V)",500,-15,5); 
    thermal_v->SetLineColor(12); 
    thermal_v->GetXaxis()->SetTitle("F"); 

     cdata_th.Draw("F>>thermal_hv", qualityCut * "isMostImpulsive * (theta<0)","goff"); 
     cdata_th.Draw("F>>thermal_h", qualityCut * "isMostImpulsive * (theta<0) * (iteration < 5)","goff"); 
     cdata_th.Draw("F>>thermal_v", qualityCut * "isMostImpulsive * (theta<0) * (iteration > 4)","goff"); 

    int nabove_hv =thermal_hv->Integral(); 
    int nabove_h = thermal_h->Integral(); 
    int nabove_v = thermal_v->Integral(); 

    TH1 * count_below = new TH1D("below","Below horizon", 2,-0.5,1.5); 
    count_below->GetYaxis()->SetTitle("vpol"); 
    cdata_th.Draw("iteration > 4 >> below", qualityCut * "isMostImpulsive * (theta > 0)","goff"); 
    int nbelow_hv = count_below->Integral(); 
    int nbelow_h = count_below->GetBinContent(1); 
    int nbelow_v = count_below->GetBinContent(2); 


    printf(" region      |  H+ V    |  H    |   V    \n"); 
    printf("-----------------------------------------\n") ; 
    printf(" above horiz | %d | %d  | %d\n", nabove_hv, nabove_h, nabove_v);  
    printf(" below horiz | %d | %d  | %d\n", nbelow_hv, nbelow_h, nbelow_v);  

    thermal_hv->Write(); 
    thermal_h->Write(); 
    thermal_v->Write(); 

    TCanvas * cthermal = new TCanvas("cthermal","Thermal",800,600); 
    thermal_hv->DrawCopy(); 
    thermal_h->DrawCopy("same"); 
    thermal_v->DrawCopy("same"); 
    cthermal->BuildLegend(); 


  }




}



void cutOptimization(bool remake_hists = true)
{

  if (remake_hists) makeHists("cutopt.root"); 
  TFile f("cutopt.root"); 

  TH2 * efficiency = (TH2*) f.Get("eff"); 
  TH1 * thermal_hv = (TH1*) f.Get("thermal_hv"); 

  double Fs[5] = {2.8,2.9,3,3.1,3.2}; 

 
  TCanvas * cthermalbg = new TCanvas("cthermalbg","Thermal BG", 1800,900); 
  TCanvas * cnobs0 = new TCanvas("cnobs0","N Observed S = 0", 1800,900); 
  TCanvas * cnobs5 = new TCanvas("cnobs5","N Observed S = 5", 1800,900); 
  TCanvas * climit = new TCanvas("climit","Limits", 1800,900); 
  TCanvas * ctotal = new TCanvas("ctotal","Total", 1800,900); 

  double average_upper_limit[5]; 
  double fraction_measured[5]; 
  double E[5]; 

  SensitivityCalculator calc; 
  calc.setAnthroBackground(1,0.97,0.19); 

  TH1 * banthro = calc.histBAnthro(); 

  TH1 * bg[5]; 
  TH1 * bg_total[5]; 
  TH1 * nobs0[5]; 
  TH1 * nobs5[5]; 


  for (int i = 0; i < 5; i++) 
  {
     
    double F = Fs[i]; 
    int bin = thermal_hv->FindBin(F); 
    int nthermal = thermal_hv->Integral(bin, thermal_hv->GetNbinsX()); 
    printf("nthermal: %g %d\n",F, nthermal); 
    calc.setThermalBackground(nthermal, 3.2); 
    E[i] = efficiency->Interpolate(12,F); 
    calc.setEfficiency( efficiency->Interpolate(12,F), 0.05); 


    bg[i] = calc.histBThermal(); 
    bg[i]->Scale(1./bg[i]->Integral()); 
    bg[i]->SetName(TString::Format("bg%d",i)); 
    bg[i]->SetTitle(TString::Format("F = %g", F)); 
    bg[i]->SetLineColor(i+1); 

    bg_total[i] = calc.histBTotal(); 
    bg_total[i]->SetName(TString::Format("bgtotal%d",i)); 
    bg_total[i]->SetTitle(TString::Format("F = %g", F)); 

    nobs0[i] = calc.histNObserved(0); 
    nobs0[i]->SetTitle(TString::Format("nobs0%d",i)); 
    nobs0[i]->SetTitle(TString::Format("F = %g", F)); 
    nobs0[i]->SetLineColor(i+1); 
    nobs5[i] = calc.histNObserved(5); 
    nobs5[i]->SetTitle(TString::Format("nobs5%d",i)); 
    nobs5[i]->SetTitle(TString::Format("F = %g", F)); 
    nobs5[i]->SetLineColor(i+1); 

    TGraphErrors * limit = calc.confidenceBands(0,20); 
    limit->SetTitle(TString::Format("F = %g\n", F)); 
    limit->SetLineColor(i+1); 
    limit->SetLineWidth(6-i); 
    limit->SetFillColor(0); 

    average_upper_limit[i] = 0; 
    fraction_measured[i] = 0; 
    for (int j = 0; j < 20; j++)
    {
      average_upper_limit[i] += (limit->GetY()[j] + limit->GetEY()[j]) * nobs0[i]->GetBinContent(j+1); 
      if ( (limit->GetY()[j] - limit->GetEY()[j])) fraction_measured[i] += nobs5[i]->GetBinContent(j+1); 
    }

    cthermalbg->cd(); 
    bg[i]->DrawCopy( i > 0 ? "same" : ""); 
    cnobs0->cd(); 
    nobs0[i]->DrawCopy( i > 0 ? "same" : ""); 
    cnobs5->cd(); 
    nobs5[i]->DrawCopy( i > 0 ? "same" : ""); 
    ctotal->cd(); 
    bg_total[i]->DrawCopy( i > 0 ? "same" : ""); 
    climit->cd(); 
    limit->Draw(i == 0 ? " ap" : "psame"); 

  }

  printf("\\begin{tabular}{l|l|l|l|l|l|l}\n"); 
  printf("F & Avg Upper Limit & S=5 Acceptance & Analysis Efficiency ($\\sigma_{\\mathcal{E}} = 0.05$) & $\\mathcal{B}_{thermal}$ & $\\mathcal{B}_{thermal+anthro}$  \\\\\n"); 

  for (int i = 0; i < 5; i++)
  {
    double qs[3] = {0.16,0.5,0.84}; 
    double qs_thermal[3], qs_total[3]; 
    bg[i]->GetQuantiles(3, qs_thermal, qs); 
    bg_total[i]->GetQuantiles(3, qs_total, qs); 
    printf("%0.2g & %0.3g & %0.3g & %0.3g & $%0.3g^{+%0.2g}_{-%0.2g}$ & $%0.3g^{+%0.2g}_{-%.02g}$ \\\\\n", Fs[i], average_upper_limit[i], fraction_measured[i], E[i] , qs_thermal[1], qs_thermal[2]-qs_thermal[1], qs_thermal[1]-qs_thermal[0], qs_total[1], qs_total[2]-qs_total[1], qs_total[1]-qs_total[0]); 
  }


}

