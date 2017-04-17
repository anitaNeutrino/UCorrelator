
#include "macro/cuts.C"
#include "sys/wait.h" 



struct pointing_info 
{

  double dphi[3][2]; // hpol/vpol/dpol  mean/rms
  double dtheta[3][2]; // hpol/vpol/dpol  mean/rms
  double fraction_good[3]; 
  double fraction_filtered[3][2]; 
}; 

struct rejection_info
{
  double wais_overlap; 
  double ldb_overlap;  
  int Nwais; 
  int Nldb; 
  int Nbg; 
}; 



  



void doPulserPlots(const char * pulser, TChain *c, pointing_info * p)
{

  const char * filter = c->GetName(); 
  gStyle->SetMarkerStyle(6); 

  ////HPOL
  TCanvas * hpol = new TCanvas(TString::Format("%s_hpol_%s", pulser,filter),TString::Format("%s H-Pol with %s",pulser,filter), 1800,900); 
  hpol->Divide(2,1); 
  hpol->cd(1); 
  TCut is_hpol_pulser("peakPulserCoherentH > 1.5 * peakPulserCoherentV && peakPulserCoherentH > 40"); 
  TCut is_vpol_pulser("peakPulserCoherentV > 1.5 * peakPulserCoherentH && peakPulserCoherentV > 40"); 
  TCut is_dpol_pulser("peakPulserCoherentV < 1.5 * peakPulserCoherentH  && peakPulserCoherentH < 1.5* peakPulserCoherentV && peakPulserCoherentV > 40 && peakPulserCoherentH >40"); 
  TCut is_bad (TString::Format("(abs(FFTtools::wrap(peak[0][0].theta-%s.theta,360,0)) > 3 || abs(FFTtools::wrap(peak[0][0].phi-%s.phi,360,0)) > 3 )" ,pulser,pulser)); 

  int Ntotal = c->Draw(TString::Format("FFTtools::wrap(peak[0][0].theta-%s.theta,360,0):FFTtools::wrap(peak[0][0].phi-%s.phi,360,0)",pulser,pulser), brightestPeak && is_hpol_pulser); 

  TString misreco_str;
  misreco_str.Form("filterPlots/misreco_%s_%s_hpol.txt",filter,pulser); 

  FILE * misreco = fopen(misreco_str.Data(),"w"); 
  int Nbad = c->Draw("run:eventNumber",is_bad && brightestPeak && is_hpol_pulser,"goff");
  for (int i = 0; i < Nbad; i++) 
  {
    int run = (int) c->GetV1()[i];
    int event = (int) c->GetV2()[i];
    fprintf(misreco,"run %d ev %d\n", run,event);
  }
  fclose(misreco); 


  hpol->cd(2); 
  TH2 * h_hpol = new TH2D(TString::Format("h_hpol_%s_%s", pulser, filter ),"Zoomed H-Pol; dphi; dtheta;",100,-3,3,100,-3,3); 

  c->Draw(TString::Format("FFTtools::wrap(peak[0][0].theta-%s.theta,360,0):FFTtools::wrap(peak[0][0].phi-%s.phi,360,0) >> h_hpol_%s_%s",pulser,pulser,pulser,filter), brightestPeak && is_hpol_pulser,"colz"); 
  int Nclose = h_hpol->Integral(); 

  p->dphi[0][0] = h_hpol->GetMean(1); 
  p->dphi[0][1] = h_hpol->GetRMS(1); 
  p->dtheta[0][0] = h_hpol->GetMean(2); 
  p->dtheta[0][1] = h_hpol->GetRMS(2); 
  p->fraction_good[0] = double(Nclose)/Ntotal; 

  hpol->SaveAs(TString::Format("filterPlots/%s_hpol_%s.png",pulser,filter)); 

  ////VPOL
  TCanvas * vpol = new TCanvas(TString::Format("%s_vpol_%s", pulser,filter),TString::Format("%s V-Pol with %s",pulser,filter), 1800,900); 
  vpol->Divide(2,1); 
  vpol->cd(1); 
  Ntotal = c->Draw(TString::Format("FFTtools::wrap(peak[1][0].theta-%s.theta,360,0):FFTtools::wrap(peak[1][0].phi-%s.phi,360,0)",pulser,pulser), brightestPeak && is_vpol_pulser); 

  misreco_str.Form("filterPlots/misreco_%s_%s_vpol.txt",filter,pulser); 
  misreco = fopen(misreco_str.Data(),"w"); 
  Nbad = c->Draw("run:eventNumber",is_bad && brightestPeak && is_vpol_pulser,"goff");
  for (int i = 0;  i < Nbad; i++) 
  {
    int run = (int) c->GetV1()[i];
    int event = (int) c->GetV2()[i];
    fprintf(misreco,"run %d ev %d\n", run,event);
  }
  fclose(misreco); 




  vpol->cd(2); 
  TH2 * h_vpol = new TH2D(TString::Format("h_vpol_%s_%s", pulser, filter ),"Zoomed V-Pol; dphi; dtheta",100,-3,3,100,-3,3); 
  Nclose = c->Draw(TString::Format("FFTtools::wrap(peak[1][0].theta-%s.theta,360,0):FFTtools::wrap(peak[1][0].phi-%s.phi,360,0) >> h_vpol_%s_%s",pulser,pulser,pulser,filter), brightestPeak && is_vpol_pulser,"colz"); 
  vpol->SaveAs(TString::Format("filterPlots/%s_vpol_%s.png",pulser,filter)); 
  Nclose = h_vpol->Integral(); 

  p->dphi[1][0] = h_vpol->GetMean(1); 
  p->dphi[1][1] = h_vpol->GetRMS(1); 
  p->dtheta[1][0] = h_vpol->GetMean(2); 
  p->dtheta[1][1] = h_vpol->GetRMS(2); 
  p->fraction_good[1] = double(Nclose)/Ntotal; 



  ////DPOL
  TCanvas * dpol = new TCanvas(TString::Format("%s_dpol_%s", pulser,filter),TString::Format("%s D-Pol with %s",pulser,filter), 1800,900); 
  dpol->Divide(2,1); 
  dpol->cd(1); 
  Ntotal = c->Draw(TString::Format("FFTtools::wrap(peak[][0].theta-%s.theta,360,0):FFTtools::wrap(peak[][0].phi-%s.phi,360,0)",pulser,pulser), brightestPeak && is_dpol_pulser); 

  misreco_str.Form("filterPlots/misreco_%s_%s_dpol.txt",filter,pulser); 
  misreco = fopen(misreco_str.Data(),"w"); 
  Nbad = c->Draw("run:eventNumber",is_bad && brightestPeak && is_dpol_pulser,"goff");
  for (int i = 0;  i < Nbad; i++) 
  {
    int run = (int) c->GetV1()[i];
    int event = (int) c->GetV2()[i];
    fprintf(misreco,"run %d ev %d\n", run,event);
  }
  fclose(misreco); 



  dpol->cd(2); 
  TH2 * h_dpol = new TH2D(TString::Format("h_dpol_%s_%s", pulser, filter ),"Zoomed D-Pol; dphi; dtheta",100,-3,3,100,-3,3); 
  Nclose = c->Draw(TString::Format("FFTtools::wrap(peak[][0].theta-%s.theta,360,0):FFTtools::wrap(peak[][0].phi-%s.phi,360,0) >> h_dpol_%s_%s",pulser,pulser,pulser,filter), brightestPeak && is_dpol_pulser,"colz"); 
  dpol->SaveAs(TString::Format("filterPlots/%s_dpol_%s.png",pulser,filter)); 

  Nclose = h_dpol->Integral(); 
  p->dphi[2][0] = h_dpol->GetMean(1); 
  p->dphi[2][1] = h_dpol->GetRMS(1); 
  p->dtheta[2][0] = h_dpol->GetMean(2); 
  p->dtheta[2][1] = h_dpol->GetRMS(2); 
  p->fraction_good[2] = double(Nclose)/Ntotal; 



  ////Filtered Fraction 
  TCanvas * fraction_filtered = new TCanvas(TString::Format("%s_%s_filter_fraction",pulser,filter), TString::Format("%s Filter Fraction with %s", pulser,filter),1800,900); 

  fraction_filtered->Divide(2,2); 
  fraction_filtered->cd(1); 

  TH2 * h_fracangle = new TH2D(TString::Format("h_fracangle_%s_%s", pulser, filter ),"Filtered Fraction by Angle (all); dphi ; dtheta",100,-10,10,100,-10,10); 
  c->Draw(TString::Format("FFTtools::wrap(peak[][].theta - %s.theta,360,0): FFTtools::wrap(peak[][].phi-%s.phi,360,0):coherent_filtered[][].totalPower / coherent[][].totalPower  >> h_fracangle_%s_%s",pulser,pulser,pulser,filter),brightestPeak,"colz"); 
  TH1 * h_frac= new TH1D(TString::Format("h_frac_%s_%s", pulser, filter ),"Filtered Fraction",100,0,1.2); 
  TH1 * h_frac_hpol= new TH1D(TString::Format("h_frac_hpol_%s_%s", pulser, filter ),"Filtered Fraction",100,0,1.2); 
  TH1 * h_frac_vpol= new TH1D(TString::Format("h_frac_vpol_%s_%s", pulser, filter ),"Filtered Fraction",100,0,1.2); 
  TH1 * h_frac_dpol= new TH1D(TString::Format("h_frac_dpol_%s_%s", pulser, filter ),"Filtered Fraction",100,0,1.2); 
  h_frac_hpol->SetLineColor(2); 
  h_frac_vpol->SetLineColor(3); 
  h_frac_dpol->SetLineColor(4); 
  fraction_filtered->cd(2); 
  c->Draw(TString::Format("coherent_filtered[][].peakHilbert / coherent[][].peakHilbert  >> h_frac_%s_%s",pulser,filter),brightestPeak,""); 


  c->Draw(TString::Format("coherent_filtered[][].peakHilbert / coherent[][].peakHilbert  >> h_frac_hpol_%s_%s",pulser,filter),brightestPeak && is_hpol_pulser,"same"); 
  p->fraction_filtered[0][0] = h_frac_hpol->GetMean(); 
  p->fraction_filtered[0][1] = h_frac_hpol->GetRMS(); 


  c->Draw(TString::Format("coherent_filtered[][].peakHilbert / coherent[][].peakHilbert  >> h_frac_vpol_%s_%s",pulser,filter),brightestPeak && is_vpol_pulser,"same"); 
  p->fraction_filtered[1][0] = h_frac_vpol->GetMean(); 
  p->fraction_filtered[1][1] = h_frac_vpol->GetRMS(); 
  c->Draw(TString::Format("coherent_filtered[][].peakHilbert / coherent[][].peakHilbert  >> h_frac_dpol_%s_%s",pulser,filter),brightestPeak && is_dpol_pulser,"same"); 

  p->fraction_filtered[2][0] = h_frac_dpol->GetMean(); 
  p->fraction_filtered[2][1] = h_frac_dpol->GetRMS(); 
  c->SetLineColor(1); 

  fraction_filtered->cd(3); 
  c->Draw(TString::Format("coherent_filtered[][].peakHilbert / coherent[][].peakHilbert: coherent.peakHilbert"),brightestPeak,"colz"); 

  fraction_filtered->cd(4); 
  c->Draw(TString::Format("coherent_filtered[][].peakHilbert / coherent[][].peakHilbert: coherent.peakHilbert"),brightestPeak,"colz"); 

  fraction_filtered->SaveAs(TString::Format("filterPlots/%s_fraction_%s.png",pulser,filter)); 

  /*
  TCanvas * distance_plot = new TCanvas(TString::Format("%s_%s_distance",filter,pulser),TString::Format("Distance Plot for %s %s",filter,pulser)); 

  TH2 * h_distance = new TH2D(TString::Format("h_distance_%s_%s", pulser, filter ),"Angular distance (all); triggerTime ; triggerTimeNs",1000,c->GetMinimum("triggerTime"),c->GetMaximum("triggerTime"),1000,0,1e9); 
  c->Draw(TString::Format("triggerTimeNs:triggerTime:sqrt(pow(FFTtools::wrap(peak[][].theta - %s.theta,360,0),2) + pow(FFTtools::wrap(peak[][].phi-%s.phi,360,0),2))  >> h_distance_%s_%s",pulser,pulser,pulser,filter),brightestPeak,"colz"); 
  distance_plot->SaveAs(TString::Format("filterPlots/%s_distance_%s.png",pulser,filter)); 
  */

  gStyle->SetMarkerStyle(1); 
}

double computeOverlap(TH2 * h, TGraph * g) 
{

  double overlap = 0; 
  for (int i =0; i < g->GetN(); i++) 
  {
    if (g->GetX()[i] < h->GetXaxis()->GetXmin() || g->GetX()[i] > h->GetXaxis()->GetXmax()) continue; 
    if (g->GetY()[i] < h->GetYaxis()->GetXmin() || g->GetY()[i] > h->GetYaxis()->GetXmax()) continue; 
    overlap += h->Interpolate(g->GetX()[i], g->GetY()[i]); 
  }
  printf("%g\n",overlap); 
  overlap /= g->GetN(); 
  overlap /= h->Integral(); 

  return overlap; 
}



void makeCutPlot(const char * filter, TChain * bg, TChain * wais, TChain * ldb, rejection_info *p) 
{

  gStyle->SetOptStat(0); 
  TCanvas * c = new TCanvas(TString::Format("ccut_%s",filter), TString::Format("Cutplot %s",filter), 1920,1080); 
  TH2I * hbg = new TH2I(TString::Format("hbg_%s", filter),TString::Format("The standard cut plot for %s; Correlation Map Peak; Coherent Peak Hilbert",filter), 100,0,0.5,100,0,300); 
  bg->Draw(TString::Format("coherent[][].peakHilbert:peak[][].value >> hbg_%s",filter),
           brightestPeak && aboveHorizon && blastCut && triggered && notMasked,"colz"); 

  p->Nbg = hbg->Integral(); 
  wais->SetMarkerColorAlpha(30,0.5); 


  int N = wais->Draw("coherent[][].peakHilbert:peak[][].value", brightestPeak && "peakPulserCoherentH > 40 && abs(FFTtools::wrap(peak[][].phi-wais.phi,360,0)) < 3 && abs(FFTtools::wrap(peak[][].theta - wais.theta,360,0)) < 3","psame"); 
  TGraph  gwais(N, wais->GetV2(), wais->GetV1()); 
  p->Nwais = N;

  ldb->SetMarkerColorAlpha(46,0.5); 
  N = ldb->Draw("coherent[][].peakHilbert:peak[][].value", brightestPeak && "(peakPulserCoherentH > 40 || peakPulserCoherentV > 40) && abs(FFTtools::wrap(peak[][].phi-ldb.phi,360,0)) < 3 && abs(FFTtools::wrap(peak[][].theta - ldb.theta,360,0)) < 3","psame"); 
  TGraph  gldb(N, ldb->GetV2(), ldb->GetV1()); 


  p->Nldb = N;
  p->ldb_overlap = computeOverlap(hbg,&gldb); 
  p->wais_overlap = computeOverlap(hbg,&gwais); 


  c->SaveAs(TString::Format("filterPlots/standard_cut_plot_%s.png",filter)); 
  gStyle->SetOptStat(1); 

}



static std::vector<int> waitforme; 
static std::vector<const char *> filters; 
static std::vector<const char *> descriptions; 

static int first = 0; 

void doFilterAlgo(const char * filter, const char * description) 
{
  pid_t pid = fork(); 

  bool output_runs = false; 
  if (!first) output_runs = true; 
  first++; 

  if (pid == 0)
  {
    TChain cwais(filter); 

    cwais.Add("filter/*_wais_*.root"); 


    pointing_info ldb_point;
    pointing_info wais_point;
    doPulserPlots("wais",&cwais,&wais_point); 

    TChain cldb(filter); 
    cldb.Add("filter/*_ldb_*.root"); 

    doPulserPlots("ldb",&cldb,&ldb_point); 

    TChain cbg(filter); 
    cbg.Add("filter/*_bg*.root"); 

    rejection_info reject; 
    makeCutPlot(filter, &cbg,&cwais,&cldb, &reject); 

    TString str(filter); 
    TString escaped = str.ReplaceAll("_","\\_"); 

    /* output the row for the WAIS summary plot */ 
    TString ofile; 
    ofile.Form("filterPlots/slides/%s_wais.row",filter); 
    FILE * row = fopen(ofile.Data(),"w"); 

    fprintf(row,"%s&%g&%g&%g&%g&%g$\\pm$ %g &%g",
                 escaped.Data(),      
                 wais_point.dphi[0][0], 
                 wais_point.dphi[0][1], 
                 wais_point.dtheta[0][0], 
                 wais_point.dtheta[0][1], 
                 wais_point.fraction_filtered[0][0], 
                 wais_point.fraction_filtered[0][1], 
                 wais_point.fraction_good[0] 
                 
        );
    fclose(row); 
                 
    /* output the row for the LDB-H summary table */ 
    ofile.Form("filterPlots/slides/%s_ldb-h.row",filter); 
    row = fopen(ofile.Data(),"w"); 

    fprintf(row,"%s&%g&%g&%g&%g&%g$\\pm$ %g &%g",
                 escaped.Data(),      
                 ldb_point.dphi[0][0], 
                 ldb_point.dphi[0][1], 
                 ldb_point.dtheta[0][0], 
                 ldb_point.dtheta[0][1], 
                 ldb_point.fraction_filtered[0][0], 
                 ldb_point.fraction_filtered[0][1], 
                 ldb_point.fraction_good[0] 
                 
        );

    fclose(row); 

   /* output the row for the LDB-V summary table */ 
    ofile.Form("filterPlots/slides/%s_ldb-v.row",filter); 
    row = fopen(ofile.Data(),"w"); 

    fprintf(row,"%s&%g&%g&%g&%g&%g$\\pm$ %g &%g",
                 escaped.Data(),      
                 ldb_point.dphi[1][0], 
                 ldb_point.dphi[1][1], 
                 ldb_point.dtheta[1][0], 
                 ldb_point.dtheta[1][1], 
                 ldb_point.fraction_filtered[1][0], 
                 ldb_point.fraction_filtered[1][1], 
                 ldb_point.fraction_good[1] 
                 
        );

    
    fclose(row); 
 
    /* output the row for the LDB-D summary table */ 
    ofile.Form("filterPlots/slides/%s_ldb-d.row",filter); 
    row = fopen(ofile.Data(),"w"); 

    fprintf(row,"%s&%g&%g&%g&%g&%g$\\pm$ %g &%g",
                   escaped.Data(),      
                   ldb_point.dphi[2][0], 
                   ldb_point.dphi[2][1], 
                   ldb_point.dtheta[2][0], 
                   ldb_point.dtheta[2][1], 
                   ldb_point.fraction_filtered[2][0], 
                   ldb_point.fraction_filtered[2][1], 
                   ldb_point.fraction_good[2] 
                   
         );

    fclose(row); 
 

    /* output the row for the background summary table */ 
    ofile.Form("filterPlots/slides/%s_background.row",filter); 
    row = fopen(ofile.Data(),"w"); 

    fprintf(row,"%s&%d&%d&%d&%g&%g",
                   escaped.Data(),      
                   reject.Nbg, reject.Nwais, reject.Nldb,
                   reject.wais_overlap, reject.ldb_overlap
         );

    fclose(row); 
 


    /* output the slides */ 

    ofile.Form("filterPlots/slides/%s.tex",filter); 
    FILE * slide = fopen(ofile.Data(),"w"); 

    fprintf(slide,"\\begin{frame}\n\\frametitle{%s - WAIS}\n\n", escaped.Data()); 
    fprintf(slide,"\\begin{columns}\n\t\\begin{column}{0.5\\textwidth}\n"); 
    fprintf(slide,"\t\t\\includegraphics[width=3in]{../wais_hpol_%s}\n", filter); 
    fprintf(slide,"\t\\end{column}\n\t\\begin{column}{0.5\\textwidth}\n"); 
    fprintf(slide,"\t\t\\includegraphics[width=3in]{../wais_fraction_%s}\n", filter); 
    fprintf(slide,"\t\\end{column}\n\\end{columns}\n\\end{frame}\n\n"); 

    fprintf(slide,"\\begin{frame}\n\\frametitle{%s - LDB Hpol / V-Pol}\n\n", escaped.Data()); 
    fprintf(slide,"\\begin{columns}\n\t\\begin{column}{0.5\\textwidth}\n"); 
    fprintf(slide,"\t\t\\includegraphics[width=3in]{../ldb_hpol_%s}\n", filter); 
    fprintf(slide,"\t\\end{column}\n\t\\begin{column}{0.5\\textwidth}\n"); 
    fprintf(slide,"\t\t\\includegraphics[width=3in]{../ldb_vpol_%s}\n", filter); 
    fprintf(slide,"\t\\end{column}\n\\end{columns}\n\\end{frame}\n\n"); 

    fprintf(slide,"\\begin{frame}\n\\frametitle{%s - LDB D-pol / Fraction}\n\n", escaped.Data()); 
    fprintf(slide,"\\begin{columns}\n\t\\begin{column}{0.5\\textwidth}\n"); 
    fprintf(slide,"\t\t\\includegraphics[width=3in]{../ldb_dpol_%s}\n", filter); 
    fprintf(slide,"\t\\end{column}\n\t\\begin{column}{0.5\\textwidth}\n"); 
    fprintf(slide,"\t\t\\includegraphics[width=3in]{../ldb_fraction_%s}\n", filter); 
    fprintf(slide,"\t\\end{column}\n\\end{columns}\n\\end{frame}\n\n"); 

    fprintf(slide,"\\begin{frame}\n\\frametitle{%s - Background Separation}\n\n", escaped.Data()); 
    fprintf(slide,"\\includegraphics[width=5in]{../standard_cut_plot_%s}\n",filter); 
    fprintf(slide,"\\end{frame}\n\n"); 

    fclose(slide); 


    if (output_runs) 
    {
      FILE * runs = fopen("filterPlots/slides/runs.tex","w"); 

      fprintf(runs,"\\begin{frame}\n\\frametitle{Runs used}\n\n\\begin{description}\n"); 
      
      int nruns = cwais.Draw("run","LocalEntry$==0","goff"); 
      fprintf(runs,"\\item[WAIS (%d runs):]", nruns); 
      for (int i = 0; i < nruns; i++) 
      {
        int run = (int) cwais.GetV1()[i]; 
        fprintf(runs,"%d ", run); 
      }
      fprintf(runs,"\n"); 

      nruns = cldb.Draw("run","LocalEntry$==0","goff"); 
      fprintf(runs,"\\item[LDB (%d runs):]", nruns); 
      for (int i = 0; i < nruns; i++) 
      {
        int run = (int) cldb.GetV1()[i]; 
        fprintf(runs,"%d ", run); 
      }
      fprintf(runs,"\n"); 

      nruns = cbg.Draw("run","LocalEntry$==0","goff"); 
      fprintf(runs,"\\item[BG (%d runs):]", nruns); 
      for (int i = 0; i < nruns; i++) 
      {
        int run = (int) cbg.GetV1()[i]; 
        fprintf(runs,"%d ", run); 
      }
      fprintf(runs,"\n"); 

      fprintf(runs,"\\end{description}\n\\end{frame}\n\n"); 

      fclose(runs); 

    }

    printf("Done with %s\n", filter); 
    _exit(0); 
  }
  else
  {
    filters.push_back(filter); 
    descriptions.push_back(description); 
    waitforme.push_back(pid); 
  }

}

void filterEvalA3Plots()
{
  system("mkdir -p filterPlots"); 
  

  doFilterAlgo("sinsub_10_0","Sine Subtract filter, 10\\% min reduction, 0 bad iters"); 
  doFilterAlgo("adsinsub_1_10_0","Adaptive Sine Subtract filter, exp=1, 10\\% min reduction, 0 bad iters"); 
  doFilterAlgo("adsinsub_2_10_0","Adaptive Sine Subtract filter, exp=2, 10\\% min reduction, 0 bad iters"); 
  doFilterAlgo("adsinsub_2_10_3","Adaptive Sine Subtract filter, exp=2, 10\\% min reduction, 3 bad iters"); 
  doFilterAlgo("adsinsub_2_10_3","Adaptive Sine Subtract filter, exp=2, 10\\% min reduction, 3 bad iters"); 
  doFilterAlgo("adsinsub_3_10_3","Adaptive Sine Subtract filter, exp=3, 10\\% min reduction, 3 bad iters"); 
  doFilterAlgo("adsinsub_2_20_0","Adaptive Sine Subtract filter, exp=2, 20\\% min reduction, 0 bad iters"); 
  doFilterAlgo("brickwall_2_0","Brickwall filter with peakiness thresh 2"); 
  doFilterAlgo("brickwall_2_1","Brickwall filter with peakiness thresh 2, filled notch"); 
  doFilterAlgo("geom","Geometric Filter (in progress)"); 

  FILE * latex = fopen("filterPlots/slides/slides.tex","w"); 
  fprintf(latex, "\\documentclass[hyperref={pdfpagelabels=false},aspectratio=169]{beamer} \\mode<presentation> { \\usetheme{Boadilla} \\usecolortheme{beaver} }\n"); 
  fprintf(latex, "\\setbeamertemplate{navigation symbols}{}\n");
  fprintf(latex, "\\title{Autogenerated Filter Plots}\n");
  time_t now = time(0); 

  char buf[80]; 
  strftime(buf,sizeof(buf),"[%m/%d/%y]{%B %d, %Y}",localtime(&now)); 
  fprintf(latex,"\\date%s\n", buf); 
  fprintf(latex,"\\author{Cosmin Deaconu}\n\n"); 
  fprintf(latex,"\\begin{document}\n"); 
  fprintf(latex,"\\begin{frame}[plain]\n\\maketitle \n\\end{frame}\n\n"); 
  fprintf(latex,"\\begin{frame}\n\\frametitle{The contenders}\n\\begin{description}\n"); 

  for (int i = 0; i < filters.size(); i++)
  {
    TString str(filters[i]); 
    TString escaped = str.ReplaceAll("_","\\_"); 
    fprintf(latex, "\t\\item[%s] %s\n", escaped.Data(), descriptions[i]); 
  }
  fprintf(latex,"\\end{description}\n\\end{frame}\n\n"); 
 
  fprintf(latex,"\\input{runs.tex}\n\n"); 

  fprintf(latex,  "\\begin{frame}\n\\frametitle{Notes}\n" 
                  "  \\begin{itemize}\n" 
                  "    \\item H-Pol means coherent peak in H-pol is 50 percent more than V-pol in pulser direction (and peak is at least 40). \n"
                  "    \\item V-Pol means coherent peak in V-pol is 50 percent more than H-pol in pulser direction (and peak is at least 40). \n"
                  "    \\item D-Pol means not H-pol or V-pol (and peak is at least 40). \n"
                  "    \\item Fraction good are events within 3 degrees in both $\\phi$ and $\\theta$\n " 
                  "    \\item Fraction filtered based on hilbert peaks. black = total, red = hpol, green = vpol, blue = dpol. \n"
                  "    \\item For background separation plot:\n"
                  "       \\begin{itemize}\n"
                  "          \\item cut on pulser events within 3 degrees of expected direction and above 40 mV hilbert peak in pulser direction.\n"
                  "          \\item  Background events from 10 percent sample where solution above zero degrees, not payload blasts, not masked, within 60 degrees of triggered phi sector \\alert{previous version of these plots avoided sun and North, not anymore}.\n" 
                  "      \\end{itemize}\n"
                  "  \\end{itemize}\n"
                  "\\end{frame}\n\n"); 

  //WAIS Summary 
  fprintf(latex,"\\begin{frame}\n\\frametitle{WAIS Pulser Summary Table}\n\\tiny\n"); 
  fprintf(latex,"\\begin{tabular}{l||c|c||c|c||c|c|}\n"); 
  fprintf(latex,"Filter & $\\bar{\\phi}$ WAIS & $\\sigma_{\\phi}$ WAIS& "
                "$ \\bar{\\theta}$ WAIS & $\\sigma_{\\theta}$ WAIS& " 
                "FilteredFrac WAIS & WAIS Purity\\\\\n\\hline\n"); 

  for (int i = 0; i < filters.size(); i++) 
  {
    fprintf(latex,"\\input{%s_wais.row}\\\\\n\\hline\n", filters[i]); 
  }

  fprintf(latex,"\\end{tabular}\n\\end{frame}\n\n"); 


  //LDB-H Summary 
  fprintf(latex,"\\begin{frame}\n\\frametitle{LDB-H Pulser Summary Table}\n\\tiny\n"); 
  fprintf(latex,"\\begin{tabular}{l||c|c||c|c||c|c|}\n"); 
  fprintf(latex,"Filter & $\\bar{\\phi}$ LDB-H & $\\sigma_{\\phi}$ LDB-H& "
                "$ \\bar{\\theta}$ LDB-H & $\\sigma_{\\theta}$ LDB-H& " 
                "FilteredFrac LDB-H & LDB-H Purity\\\\\n\\hline\n"); 

  for (int i = 0; i < filters.size(); i++) 
  {
    fprintf(latex, "\\input{%s_ldb-h.row}\\\\\n\\hline\n", filters[i]); 
  }

  fprintf(latex,"\\end{tabular}\n\\end{frame}\n\n"); 

  //LDB-V Summary 
  fprintf(latex,"\\begin{frame}\n\\frametitle{LDB-V Pulser Summary Table}\n\\tiny\n"); 
  fprintf(latex,"\\begin{tabular}{l||c|c||c|c||c|c|}\n"); 
  fprintf(latex,"Filter & $\\bar{\\phi}$ LDB-V & $\\sigma_{\\phi}$ LDB-V& "
                "$ \\bar{\\theta}$ LDB-V & $\\sigma_{\\theta}$ LDB-V& " 
                "FilteredFrac LDB-V & LDB-V Purity\\\\\n\\hline\n"); 

  for (int i = 0; i < filters.size(); i++) 
  {
    fprintf(latex,"\\input{%s_ldb-v.row}\\\\\n\\hline\n", filters[i]); 
  }

  fprintf(latex,"\\end{tabular}\n\\end{frame}\n\n"); 

  //LDB-D Summary 
  fprintf(latex,"\\begin{frame}\n\\frametitle{LDB-D Pulser Summary Table}\n\\tiny\n"); 
  fprintf(latex,"\\begin{tabular}{l||c|c||c|c||c|c|}\n"); 
  fprintf(latex,"Filter & $\\bar{\\phi}$ LDB-D & $\\sigma_{\\phi}$ LDB-D& "
                "$ \\bar{\\theta}$ LDB-D & $\\sigma_{\\theta}$ LDB-D& " 
                "FilteredFrac LDB-D & LDB-D Purity\\\\\n\\hline\n"); 

  for (int i = 0; i < filters.size(); i++) 
  {
    fprintf(latex,"\\input{%s_ldb-d.row}\\\\\n\\hline\n", filters[i]); 
  }

  fprintf(latex,"\\end{tabular}\n\\end{frame}\n\n"); 

  //Background Separation Summary Table
  fprintf(latex,"\\begin{frame}\n\\frametitle{Background Separation Summary Table}\n\\tiny\n"); 
  fprintf(latex,"\\begin{tabular}{l||c|c|c||c|c|}\n"); 
  fprintf(latex,"Filter & NBg & NWais & NLDB & WAIS Overlap & LDB Overlap \\\\\n\\hline\n"); 
  for (int i = 0; i < filters.size(); i++) 
  {
    fprintf(latex, "\\input{%s_background.row}\\\\\n\\hline\n", filters[i]); 
  }

  fprintf(latex,"\\end{tabular}\n\\end{frame}\n\n"); 



  for (int i = 0 ; i < filters.size(); i++) 
  {
    fprintf(latex,"\\input{%s.tex}\n\n", filters[i]); 
  }

  fprintf(latex,"\\end{document}\n"); 
  fclose(latex); 
  


  //reap
  for (unsigned i = 0; i < waitforme.size(); i++)
  {
    int dummy; 
    waitpid(waitforme[i],&dummy,0); 
  }

  //make latex! 

  chdir("filterPlots/slides"); 
  system("pdflatex slides.tex"); 
  system("pdflatex slides.tex"); 

}


