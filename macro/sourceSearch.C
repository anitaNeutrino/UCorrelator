#include "UCKDE.h" 



double getDeltaPhi(double pol) 
{

  static TMVA::Reader * reader = 0; 

  static float polAngle, polAngle2, polAngle3, polAngle4, polAngle5; 
  static float blahF; 
  static int blahI; 

  if (!reader) 
  {
    reader = new TMVA::Reader(); 
    reader->AddVariable("polAngle", &polAngle);
    reader->AddVariable("polAngle2",&polAngle2); 
    reader->AddVariable("polAngle3",&polAngle3); 
    reader->AddVariable("polAngle4",&polAngle4); 
    reader->AddVariable("polAngle5",&polAngle5); 
    reader->AddSpectator("pointsToMC",&blahI); 
    reader->AddSpectator("isMostImpulsive",&blahI); 
    reader->AddSpectator("weight",&blahF); 
    reader->BookMVA("ldphi","dphi_training/weights/pointing_LD.weights.xml"); 
  }

  polAngle = pol; 
  polAngle2= pol*pol; 
  polAngle3= polAngle2*pol; 
  polAngle4= polAngle2*polAngle2; 
  polAngle5= polAngle4*pol; 
   
  return reader->EvaluateMVA("ldphi"); 
}


double getDeltaTheta(double pol, double the_theta, double lowness, double slope) 
{
  static TMVA::Reader * reader = 0; 
  static float polAngle, polAngle2, polAngle3, polAngle4, polAngle5; 
  static float theta;
  static float bandwidthMeasure, bandwidthMeasure2, bandwidthMeasure3;
  static float spectralSlope, spectralSlope2, spectralSlope3;
  static int blahI; 
  static float blahF;

  if (!reader) 
  {
    reader = new TMVA::Reader(); 
    reader->AddVariable("theta",&theta); 
    reader->AddVariable("polAngle", &polAngle);
    reader->AddVariable("polAngle2",&polAngle2); 
    reader->AddVariable("polAngle3",&polAngle3); 
    reader->AddVariable("polAngle4",&polAngle4); 
    reader->AddVariable("polAngle5",&polAngle5); 
    reader->AddVariable("bandwidthMeasure", &bandwidthMeasure);
    reader->AddVariable("bandwidthMeasure2",&bandwidthMeasure2); 
    reader->AddVariable("bandwidthMeasure3",&bandwidthMeasure3); 
    reader->AddVariable("spectralSlope", &spectralSlope);
    reader->AddVariable("spectralSlope2",&spectralSlope2);
    reader->AddVariable("spectralSlope3",&spectralSlope3); 
    reader->AddSpectator("pointsToMC",&blahI); 
    reader->AddSpectator("isMostImpulsive",&blahI); 
    reader->AddSpectator("weight",&blahF); 

    reader->BookMVA("ld","dtheta_training/weights/pointing_LD.weights.xml"); 
  }

  polAngle = pol; 
  polAngle2= pol*pol; 
  polAngle3= polAngle2*pol; 
  polAngle4= polAngle2*polAngle2; 
  polAngle5= polAngle4*pol; 
  theta = the_theta; 
  bandwidthMeasure = lowness; 
  bandwidthMeasure2 = lowness*lowness; 
  bandwidthMeasure3 = lowness*lowness*lowness; 
  spectralSlope = slope; 
  spectralSlope2 = slope*slope; 
  spectralSlope3 = slope*slope*slope; 
 
  return reader->EvaluateMVA("ld") ;
}


void getRADec(double & ra, double &dec, const AnitaEventSummary * sum, const UsefulAdu5Pat * pat,  int pol = 1, int i = 0, TH1 * dtheta_hist = 0, TH1 * dphi_hist = 0)
{
  double theta = sum->peak[pol][i].theta;
  double phi = sum->peak[pol][i].phi;
  double polang = FFTtools::wrap(90/TMath::Pi()*TMath::ATan2(sum->deconvolved[pol][i].U, sum->deconvolved[pol][i].Q),180); 
  double lowness = sum->deconvolved[pol][i].bandwidthMeasure; 
  double spectralSlope = sum->deconvolved[pol][i].spectrumSlope; 
  double dtheta = getDeltaTheta(polang,theta,lowness,spectralSlope); 
  double dphi = getDeltaPhi(polang); 
  if (dtheta_hist) 
  {
    dtheta_hist->Fill( -(dtheta-sum->peak[pol][i].theta) - sum->mc.nuTheta, sum->mc.weight); 
  }
  if (dphi_hist)
  {
    dphi_hist->Fill( FFTtools::wrap( -(dphi-sum->peak[pol][i].phi) - sum->mc.nuPhi,360,0), sum->mc.weight); 
  }
  pat->astronomicalCoordinates(phi - dphi, theta-dtheta, &ra, &dec); 
  ra = FFTtools::wrap(ra,24); 
}



static const char * the_source_distances_file = "source_distances.root"; 
static const char * the_data_overlaps_file = "A3overlaps.root"; 



double doInterpolate(TH3 * h, double RA, double dec, double unixTime, bool extend = false)
{

  if (!h) return -1;

  if ( h->GetZaxis()->GetXmin() > unixTime || h->GetZaxis()->GetXmax() < unixTime) 
  {
    if (extend) 
    {

      if (h->GetZaxis()->GetXmin() > unixTime) unixTime = h->GetZaxis()->GetXmin()+1;
      else unixTime = h->GetXaxis()->GetXmax()-1; 

    }
    else
    {
    // printf("OOB: %f %f %f %f %f\n", h->GetZaxis()->GetXmin(), h->GetZaxis()->GetXmax(), unixTime, RA, dec); 
     return 1; 
    }
  }

  int bin =  h->FindFixBin(RA,dec, unixTime); 
  return h->GetBinContent(bin); 
}


void fillOverlapHists(const char * addstring, const char * outfile, bool mc, int max = -1) 
{
  TChain c("overlap"); 
  c.Add(addstring); 

  double O; 
  double F; 
  int run; 
  int event; 
  int pol; 
  int peak; 
  double weight; 


  c.SetBranchAddress("O",&O); 
  c.SetBranchAddress("F",&F); 
  c.SetBranchAddress("run",&run); 
  c.SetBranchAddress("event",&event); 
  c.SetBranchAddress("pol",&pol); 
  c.SetBranchAddress("peak",&peak); 
  c.SetBranchAddress("weight",&weight); 

  AnitaTemplateSummary ats; 
  AnitaTemplateMachine *atm = new AnitaTemplateMachine; 
  atm->loadTemplates(); 


  TFile f(outfile,"RECREATE");  
  TTree *  out = new TTree("overlap","overlap"); 
  double RA, dec, unixTime; 
  double D_grb=0, D_txs=0, D_fava=0, D_sn=0;
  double RA_CR, dec_CR; 
  double cray_4; 
  out->Branch("O",&O); 
  out->Branch("F",&F); 
  out->Branch("run",&run); 
  out->Branch("event",&event); 
  out->Branch("pol",&pol); 
  out->Branch("peak",&peak); 
  out->Branch("RA",&RA); 
  out->Branch("dec",&dec); 

  out->Branch("RA_CR",&RA_CR); 
  out->Branch("dec_CR",&dec_CR); 
  out->Branch("unixTime",&unixTime); 
  out->Branch("weight",&weight); 
  out->Branch("cray_4",&cray_4); 

  out->Branch("D_grb",&D_grb); 
  out->Branch("D_txs",&D_txs); 
  out->Branch("D_fava",&D_fava); 
  out->Branch("D_sn",&D_sn); 

  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 3; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
  cfg.delay_to_center = true; 
  cfg.r_time_shift_correction = true; 

  cfg.combine_nantennas = 15; 
  cfg.zoomed_nant = 15; 
  cfg.use_coherent_spectra = true; 
  cfg.cross_correlate_hv = 1; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 
  FilterStrategy * strat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2"); 

  UCorrelator::setAdaptiveFilterSpectrumAverageNSecs(10); 


  //try to load the histograms 

  TFile fsd(the_source_distances_file); 

  TH3 * h_grb = (TH3*) fsd.Get("grb"); 
  TH3 * h_txs = (TH3*) fsd.Get("txs"); 
  TH3 * h_fava = (TH3*) fsd.Get("fava_equal"); 
  TH3 * h_sn = (TH3*) fsd.Get("sn"); 


  int last_run =-1; 
  AnitaDataset *d = 0;

  for (int i = 0; i < c.GetEntries(); i++) 
  {
    if (max > 0 && i > max) break; 
    c.GetEntry(i); 
    if (O <= -1)
    {
      if (!mc) continue; 
      O=-1; 
    }
    if (O > 1 && !mc) continue; 
    if (O < 0 && O > -1) O = 0; 
    if (!d || run !=last_run) 
    {
      if (d) delete d; 
      d = new AnitaDataset (run,false,WaveCalType::kDefault,mc ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA3_ROOT_DATA); 
      last_run = run; 
    }
    unixTime = d->header()->triggerTime + 1e-9 * d->header()->triggerTimeNs;
    d->getEvent(event); 
    FilteredAnitaEvent ev(d->useful(),strat, d->gps(), d->header()); 
    AnitaEventSummary sum; 
    analyzer->analyze(&ev,&sum,d->truth()); 

    atm->doTemplateAnalysis( analyzer->getCoherent(AnitaPol::AnitaPol_t(pol),peak),pol,peak,&ats);
    cray_4 = ats.coherent[pol][peak].cRay[4]; 

    UsefulAdu5Pat gps(d->gps()); 
    gps.astronomicalCoordinates(sum.peak[pol][peak].phi, sum.peak[pol][peak].theta, &RA_CR, &dec_CR); 

    getRADec(RA,dec, &sum, &gps, pol, peak); 


    D_grb =doInterpolate(h_grb, RA,dec,unixTime);
    D_txs = doInterpolate(h_txs,RA,dec,unixTime,true); 
    D_fava = doInterpolate(h_fava,RA,dec,unixTime,true); 
    D_sn = doInterpolate(h_sn,RA,dec,unixTime,true); 



    f.cd(); 
    out->Fill(); 
  }
  f.cd(); 
  out->Write(); 
} 


const int source_n_ra_bins = 720;
const int source_n_dec_bins = 360;





UCorrelator::KDE2D * fillTMVARADec(TChain * c, double mintime =0, double maxtime = 0, TH2 * h = 0, TH1* dphi_hist = 0, TH1 * dtheta_hist = 0, 
                                   const char * summary_name = "summary", const char * pat_name = "pat", 
                                   bool mostImpulsive = true, double mc_phi_cut = 5.5, double mc_theta_cut = 3.5)
{

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * pat = new Adu5Pat; 

  c->SetBranchAddress(summary_name, &sum); 
  c->SetBranchAddress(pat_name, &pat); 

  std::vector<double> ras; 
  std::vector<double> decs; 
  std::vector<double> ws; 
  printf("min_time: %g, max_time: %g, delta: %g\n", mintime, maxtime, maxtime-mintime); 
  for (int i = 0; i < c->GetEntries(); i++) 
  {
    c->GetEntry(i); 
    if (mintime && pat->realTime < mintime) continue; 
    if (maxtime && pat->realTime > maxtime) continue; 
    int pol = sum->mostImpulsivePolAsInt(); 
    int idx = sum->mostImpulsiveInd(); 

    if (abs(FFTtools::wrap(sum->peak[pol][idx].theta-sum->mc.theta,360,0)) < mc_theta_cut && abs(FFTtools::wrap(sum->peak[pol][idx].phi-sum->mc.phi,360,0)) < mc_phi_cut)
    {
       double ra,dec; 
       UsefulAdu5Pat gps(pat); 
       getRADec(ra,dec, sum, &gps, pol, idx, dtheta_hist, dphi_hist); 
       printf("%d %g %g %g %g\n",i, ra,dec, sum->peak[pol][idx].phi, sum->peak[pol][idx].theta); 
       if (h) h->Fill(ra,dec, sum->mc.weight); 
       ras.push_back(ra); 
       decs.push_back(dec); 
       ws.push_back(sum->mc.weight); 
    }
  }


  delete sum; 
  delete pat; 
  UCorrelator::KDE2D::KDE2DOptions opt(0,0.5,0.5); 
  return new UCorrelator::KDE2D(ras.size(), &ras[0], &decs[0], &ws[0], opt); 
}





void makeSourceDistanceHists(const char * dir, const char * key, const char * time_file, int start_run, int end_run) 
{

  TChain c("simulation"); 
  for (int i = start_run; i < end_run; i++) c.Add(TString::Format("%s/%d_*.root",dir, i)); 

  double start = c.GetMinimum("realTime");
  double end = c.GetMaximum("realTime");

  //figure out which times are relevant

  std::vector<double> time_bins; 
  time_bins.push_back(start); 

  printf("%f\n",start);
  if (time_file) 
  {
    FILE * tf = fopen(time_file,"r"); 
    while (!feof(tf))
    {
      double t; 
      fscanf(tf,"%lf\n", &t); 
      printf("%f\n",t); 
      if (t <= start)  continue; 
      if (t >= end) continue; 
      time_bins.push_back(t); 
    }
  }

  printf("%f\n",end);
  time_bins.push_back(end); 

  std::vector<double> ra_bins(source_n_ra_bins+1); 
  std::vector<double> dec_bins(source_n_dec_bins+1); 

  double d_ra = 24./source_n_ra_bins;
  double d_dec = 180./source_n_dec_bins;

  for (int i = 0; i <= source_n_ra_bins; i++) ra_bins[i] = d_ra*i;
  for (int i = 0; i <= source_n_dec_bins; i++) dec_bins[i] = d_dec*i-90;


  TFile f(the_source_distances_file,"UPDATE"); 
  TH3 * h = new TH3F(key,key, source_n_ra_bins, &ra_bins[0], source_n_dec_bins, &dec_bins[0], time_bins.size()-1, &time_bins[0]); 

  //ok now get the kde's 

  h->SetEntries(h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ()); 
  for (unsigned ti = 1; ti < time_bins.size(); ti++) 
  {
    TH2D this_hist(TString::Format("hradec_%d",ti), TString::Format("RA dec %d", ti), 360, 0 , 24, 180,-90,90); 
    UCorrelator::KDE2D * kde = fillTMVARADec(&c, time_bins[ti-1], time_bins[ti], &this_hist); 

    if (kde->getN() == 0) 
    {
      printf("No events! Filling histogram with 1"); 

      for (int yi = 1; yi <= h->GetNbinsY(); yi++) 
      {
        for (int xi = 1; xi<= h->GetNbinsX(); xi++) 
        {
          h->SetBinContent(xi,yi,ti, 1);
        }
      }
      delete kde;
      continue;
    }

    //now let's fill the  histogram! 

    TH2D htmp("htmp", "htmp", source_n_ra_bins, &ra_bins[0], source_n_dec_bins, &dec_bins[0]);

    for (int yi = 1; yi <= h->GetNbinsY(); yi++) 
    {
      for (int xi = 1; xi <= h->GetNbinsX(); xi++) 
      {
        double ra = h->GetXaxis()->GetBinCenter(xi); 
        double dec = h->GetYaxis()->GetBinCenter(yi); 
       //get wrap around right! 
        htmp.SetBinContent(xi,yi, kde->operator()(ra,dec) + kde->operator()(ra+24,dec) + kde->operator()(ra-24,dec));
      }
    }

    delete kde; 

    TCanvas * can = new TCanvas; 
    can->Divide(3,1);
    can->cd(1);
    htmp.DrawCopy("col2z"); 
    TH2 * pct = UCorrelator::image::makePctileHist(&htmp); 
    can->cd(2);
    pct->DrawCopy("col2z"); 
    for (int yi = 1; yi <= h->GetNbinsY(); yi++) 
    {
      for (int xi = 1; xi<= h->GetNbinsX(); xi++) 
      {
        h->SetBinContent(xi,yi,ti, pct->GetBinContent(xi,yi));
      }
    }

    can->cd(3); 
    this_hist.DrawCopy("col2z"); 
    delete pct; 
  }

  f.cd();
  h->Write(h->GetName(), TObject::kOverwrite); 
}

const double D_max = 0.975; 


TGraph * weightedAreas(TH3 * source_distances, TH1 * dec_weight = 0) 
{

  TGraph * g = new TGraph; 
  g->SetTitle("Area Ratios"); 

  double D = 0.025;
  while (D <= 1.01 ) 
  {
    g->SetPoint(g->GetN(), D, 0); 
    D+=0.025; 
  }


  for (int t = 1; t <= source_distances->GetNbinsZ(); t++) 
  {

    double dt = source_distances->GetZaxis()->GetBinWidth(t) / ( source_distances->GetZaxis()->GetXmax() - source_distances->GetZaxis()->GetXmin()); 

    for (int d = 0; d < g->GetN(); d++) 
    {
      D = g->GetX()[d]; 
      for (int j = 1; j < source_distances->GetNbinsY(); j++) 
      {
        double weight = dec_weight ? dec_weight->Interpolate(source_distances->GetYaxis()->GetBinCenter(j)) : 1; 
        for (int i = 1; i < source_distances->GetNbinsX(); i++) 
        {
          if ( source_distances->GetBinContent(i,j,t) <= D ) g->GetY()[d] += dt * weight;  
        }
      }
    }
  }

  double max = g->GetY()[g->GetN()-1]; 
  for (int i = 0; i < g->GetN(); i++)g->GetY()[i]/=max; 
  return g; 
}

double getUpperRightSum(TH2 * h, int start_i, int start_j, bool include_overlow) 
{

  double sum = 0; 
  for (int j = start_j; j <=h->GetNbinsY() +(include_overlow ? 1 : 0); j++) 
  {
    for (int i = start_i; i <=h->GetNbinsX() +(include_overlow ? 1 : 0); i++) 
    {
        sum += h->GetBinContent(i,j); 
    }
  }


       
  return sum; 

}

TH3D * getEfficiency(TTree * overlap_tree, const char * thermal_tree_file,  const char * key, double fraction = 1) 
{


  
   
  //get the denominator from the thermal tree
 

  TFile fth(thermal_tree_file); 
  TTree * th = (TTree*) fth.Get("simulation"); 

  th->Draw("F","weight * isMostImpulsive","goff"); 

  TH1 * htemp = (TH1*) gDirectory->Get("htemp");
  double denom = htemp->Integral()*fraction; 
  printf("Denom: %g\n", denom); 

  fth.Close(); 
  gROOT->cd(); 

  int N = overlap_tree->Draw(TString::Format("F:O:D_%s",key),"weight","goff"); 

  double *W = overlap_tree->GetW(); 
  double *vF = overlap_tree->GetV1();
  double *vO = overlap_tree->GetV2();
  double *vD = overlap_tree->GetV3();

  TH3D * eff = new TH3D(TString::Format("%s__eff", key), TString::Format("Efficiency Estimate for %s", key), 
                       30, 0, 3,  //F
                       120, 0, 12,  //log O 
                       40, 0, 1);  // D




  for (int k =1;  k <= eff->GetNbinsZ(); k++) 
  {
    double D =eff->GetZaxis()->GetBinLowEdge(k);
    for (int j = 1; j <= eff->GetNbinsY(); j++) 
    {
      double O = exp(-eff->GetYaxis()->GetBinLowEdge(j)); 
      for (int i = 1; i <= eff->GetNbinsX(); i++) 
      {
        double F = eff->GetXaxis()->GetBinLowEdge(i); 

        double num = 0;
        for (int ii = 0; ii < N; ii++) 
        {
          if (vD[ii] > D) continue; 
          if (vO[ii] > O) continue; 
          if (vO[ii] <=-1) continue; 
          if (vF[ii] < F) continue; 
          num += W[ii]; 
        }
        printf("%g %g %g %g\n", D, -log(O), F, num); 

        eff->SetBinContent(i,j,k,num/denom); 
      }
    }
  }

  //delete tmp; 
  return eff; 

}

TH3D * getBackgroundEstimate(TTree * overlap_tree, TH3 * source_distances, const char * key) 
{


  //get the weight histogram 
  
  TH1 * weights = new TH1D("bg_weights","Background weights",  90,-90,90); 
  overlap_tree->Draw("dec >> bg_weights","pol==1"); 

 
  
  // calculate areas of each source distacnce
  TGraph * areas = weightedAreas(source_distances, weights); 
  
  areas->Draw("alp"); 
 

  TH2 * tmp = new TH2D(TString::Format("%s_tmp_bg", key), TString::Format("Differential BG Estimate for %s", key), 
                       30, 0, 3,  //F
                       120, 0, 12);  //log O 


  //this is the differential background estimate, without the scaling
  overlap_tree->Draw(TString::Format("-log(O+1e-20):F >> %s_tmp_bg", key), TString::Format("pol==1 && D_%s >= 0.95", key),"goff"); 
  
   
  TH3D * bg = new TH3D(TString::Format("%s__bg", key), TString::Format("BG Estimate for %s", key), 
                       30, 0, 3,  //F
                       120, 0, 12,  //log O 
                       40, 0, 1);  // D



  //now we need to scale by the fraction of areas 

  for (int k =1;  k <= bg->GetNbinsZ(); k++) 
  {
    double scale = areas->Eval(bg->GetZaxis()->GetBinLowEdge(k)) / (1-areas->Eval(D_max)); 
    for (int j = 1; j <= bg->GetNbinsY(); j++) 
    {
      for (int i = 1; i <= bg->GetNbinsX(); i++) 
      {

        double sum = getUpperRightSum(tmp,i,j,true); 
        printf(" D=%g -log O =%g  F=%g sum = %g  scale =%g\n",bg->GetZaxis()->GetBinLowEdge(k), bg->GetYaxis()->GetBinLowEdge(j), bg->GetXaxis()->GetBinLowEdge(i), sum,scale); 
        bg->SetBinContent(i,j,k, sum*scale); //jefrey's prior 
        bg->SetBinError(i,j,k, (0.5+sqrt(sum+0.25))*scale);
      }
    }
  }

  //delete tmp; 
  return bg; 

}

