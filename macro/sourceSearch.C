#include "UCKDE.h" 
#include "AnitaTMVA.h" 
#include "TH3.h" 
#include "TFile.h" 
#include "TGraph.h" 
#include "AnitaEventSummary.h" 
#include "UsefulAdu5Pat.h" 
#include "TChain.h" 
#include "TFeldmanCousins.h" 
#include "AnitaTemplates.h" 
#include "AnalysisConfig.h" 
#include "Analyzer.h" 
#include "UCImageTools.h" 
#include "UCKDE.h" 
#include "FilterStrategy.h" 
#include "FilteredAnitaEvent.h" 
#include "SineSubtract.h" 
#include "AnitaDataset.h" 
#include "UCFilters.h" 


const int noffsets = 100; 
double offset_min = 3600*1.5; 
double offset_max = 3600*22.5; 

double getDeltaPhi(double pol, double the_theta) 
{

  static TMVA::Reader * reader = 0; 

  static float polAngle, polAngle2, polAngle3, polAngle4, polAngle5; 
  static float theta; 
  static float blahF; 
  static int blahI; 

  if (!reader) 
  {
    reader = new TMVA::Reader(); 
    reader->AddSpectator("dtheta",&blahF); 
    reader->AddSpectator("pointsToMC",&blahI); 
    reader->AddSpectator("isMostImpulsive",&blahI); 
    reader->AddVariable("theta",&theta); 
    reader->AddSpectator("weight",&blahF); 
    reader->AddSpectator("energy",&blahF); 
    reader->AddSpectator("peakHilbert",&blahF); 
    reader->AddSpectator("phi",&blahF); 
    reader->AddSpectator("run",&blahI); 
    reader->AddSpectator("localentry",&blahI); 
    reader->AddVariable("polAngle", &polAngle);
    reader->AddVariable("polAngle2",&polAngle2); 
    reader->AddVariable("polAngle3",&polAngle3); 
    reader->AddVariable("polAngle4",&polAngle4); 
    reader->AddVariable("polAngle5",&polAngle5); 
 
    reader->BookMVA("ldphi","dphi_training/weights/pointing_LD.weights.xml"); 
  }

  polAngle = pol; 
  polAngle2= pol*pol; 
  polAngle3= polAngle2*pol; 
  polAngle4= polAngle2*polAngle2; 
  polAngle5= polAngle4*pol; 
  theta = the_theta; 
   
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
    reader->AddSpectator("energy",&blahF); 
    reader->AddSpectator("peakHilbert",&blahF); 
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
  double dphi = getDeltaPhi(polang,theta); 
  if (dtheta_hist) 
  {
    dtheta_hist->Fill( -(dtheta-sum->peak[pol][i].theta) - sum->mc.nuTheta, sum->mc.weight); 
  }
  if (dphi_hist)
  {
    dphi_hist->Fill( FFTtools::wrap( -(dphi-sum->peak[pol][i].phi) - sum->mc.nuPhi,360,0), sum->mc.weight); 
  }
  pat->astronomicalCoordinates(phi - dphi, theta-dtheta, &ra, &dec); 
  if (ra > 24) ra-=24; 
  if (ra < 0) ra+=24; 
}



static const char * the_source_distances_file = "source_distances.root"; 
static const char * the_data_overlaps_file = "A3overlaps.root"; 



double doInterpolate(TH3 * h, double RA, double dec, double unixTime, bool extend = false)
{

  if (!h) return -1;
  if (RA > 24) RA-=24; 
  if (RA < 0) RA+=24; 

  if ( h->GetZaxis()->GetXmin() > unixTime || h->GetZaxis()->GetXmax() < unixTime) 
  {
    if (extend) 
    {

      if (h->GetZaxis()->GetXmin() > unixTime) unixTime = h->GetZaxis()->GetXmin()+1;
      else unixTime = h->GetZaxis()->GetXmax()-1; 

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


void fillOverlapHists(const char * addstring, const char * outfile, const char * mc,  int max = -1) 
{
  TChain c("overlap"); 
  c.Add(addstring); 
  printf("nentries: %lld\n", c.GetEntries()); 

  double O; 
  double F; 
  int run; 
  int ievent;
  int buffer;
  unsigned event; 
  int buffer2; 
  int pol; 
  int peak; 
  double weight = 0; 
  double mcE = 0; 


  c.SetBranchAddress("O",&O); 
  c.SetBranchAddress("F",&F); 
  c.SetBranchAddress("run",&run); 
  c.SetBranchAddress("mcE",&mcE); 
  if (mc) 
  {
    c.SetBranchAddress("event",&event); 
  }
  else
  {
    c.SetBranchAddress("event",&ievent); 
  }
  c.SetBranchAddress("pol",&pol); 
  c.SetBranchAddress("peak",&peak); 
  c.SetBranchAddress("weight",&weight); 
  c.GetEntry(0); 

  AnitaTemplateSummary ats; 
  AnitaTemplateMachine *atm = new AnitaTemplateMachine; 
  if (!mc) 
    atm->loadTemplates(); 


  TFile f(outfile,"RECREATE");  
  TTree *  out = new TTree("overlap","overlap"); 
  double RA, dec, unixTime; 


  double D_grb=0, D_txs=0, D_fava=0, D_sn=0, D_grb_12 = 0, D_grb_24 = 0;
  double D_grb_offset[noffsets] = {0}; 
  double D_grb_12_offset[noffsets] = {0}; 
  double D_grb_24_offset[noffsets] = {0}; 
  double D_txs_offset[noffsets] = {0}; 
  double D_fava_offset[noffsets] = {0}; 
  double D_sn_offset[noffsets] = {0}; 
  double offsets[noffsets] = {0}; 

  double cohSNR = 0; 
  double RA_CR = 0, dec_CR = 0; 
  double cray_4 = 0; 

  out->Branch("O",&O); 
  out->Branch("F",&F); 
  out->Branch("run",&run); 
  out->Branch("event",&event,"event/i"); 
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
  out->Branch("D_grb_12",&D_grb_12); 
  out->Branch("D_grb_24",&D_grb_24); 
  out->Branch("D_txs",&D_txs); 
  out->Branch("D_fava",&D_fava); 
  out->Branch("D_sn",&D_sn); 
  out->Branch("mcE",&mcE); 
  out->Branch("cohSNR",&cohSNR); 

  int noff = noffsets; 
  out->Branch("noffsets",&noff); 

  out->Branch("D_grb_offset",D_grb_offset,"D_grb_offset[noffsets]/D"); 
  out->Branch("D_grb_12_offset",D_grb_12_offset,"D_grb_12_offset[noffsets]/D"); 
  out->Branch("D_grb_24_offset",D_grb_24_offset,"D_grb_24_offset[noffsets]/D"); 
  out->Branch("D_txs_offset",D_txs_offset,"D_txs_offset[noffsets]/D"); 
  out->Branch("D_fava_offset",D_fava_offset,"D_fava_offset[noffsets]/D"); 
  out->Branch("D_sn_offset",D_sn_offset,"D_sn_offset[noffsets]/D"); 
  out->Branch("offsets",offsets,"offsets[noffsets]/D"); 



  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 3; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
  cfg.delay_to_center = true; 
  cfg.r_time_shift_correction = true; 

  cfg.combine_nantennas = 15; 
  cfg.zoomed_nant = 15; 
  cfg.use_coherent_spectra = true; 
  cfg.cross_correlate_hv = 0; 

  UCorrelator::Analyzer * analyzer = mc ? 0 :  new UCorrelator::Analyzer(&cfg,true) ;
  FilterStrategy * strat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2"); 

  UCorrelator::setAdaptiveFilterSpectrumAverageNSecs(10); 

  // sumfile, sumtree for mc (otherwise, we recalculate) 

  TFile * sumfile = 0;
  TTree * sumtree = 0;

  //try to load the histograms 

  TFile fsd(the_source_distances_file); 

  TH3 * h_grb = (TH3*) fsd.Get("grb"); 
  TH3 * h_grb_12 = (TH3*) fsd.Get("grb_12"); 
  TH3 * h_grb_24 = (TH3*) fsd.Get("grb_24"); 
  TH3 * h_txs = (TH3*) fsd.Get("txs"); 
  TH3 * h_fava = (TH3*) fsd.Get("fava"); 
  TH3 * h_sn = (TH3*) fsd.Get("sn"); 


  int last_run =-1; 
  AnitaDataset *d = 0;

  AnitaEventSummary sum; 
  AnitaEventSummary *psum = &sum; 
  for (int i = 0; i < c.GetEntries(); i++) 
  {
    if (max > 0 && i > max) break; 
    c.GetEntry(i); 
    if (!mc) event = ievent; 
    if (O <= -1)
    {
      if (!mc) continue; 
      O=-1; 
    }
    if (O > 10) 
    {
      if (!mc && gRandom->Uniform(0,1) > 10./O) continue; //only keep a small fraction of these 
    }

    if (O < 0 && O > -1) O = 0; 
    if (!d || run !=last_run) 
    {
      printf("run: %d\n", run); 
      if (d) delete d; 
      d = new AnitaDataset (run,false,WaveCalType::kDefault,mc ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA3_ROOT_DATA); 
      last_run = run; 
      if (mc) 
      {
        if (sumfile) delete sumfile;
        sumfile = new TFile(Form("%s/%d_sinsub_10_3_ad_2.root", mc,run)); 
        sumtree= (TTree*) sumfile->Get("simulation"); 
        sumtree->SetBranchAddress("summary",&psum); 
        sumtree->BuildIndex("run","eventNumber % (1 << 31)  "); 

      }
    }
    d->getEvent(event); 
    unixTime = d->header()->triggerTime + 1e-9 * d->header()->triggerTimeNs;
    if (!mc) 
    {
      FilteredAnitaEvent ev(d->useful(),strat, d->gps(), d->header()); 
      analyzer->analyze(&ev,&sum,d->truth()); 

      atm->doTemplateAnalysis( analyzer->getCoherent(AnitaPol::AnitaPol_t(pol),peak),pol,peak,&ats);
      cray_4 = ats.coherent[pol][peak].cRay[4]; 
    }
    else
    {
      sumtree->GetEntryWithIndex(run,event % (1<<31)); 
    }
    printf("%u %u %u\n", event, d->header()->eventNumber, sum.eventNumber); 

    UsefulAdu5Pat gps(d->gps()); 
    gps.astronomicalCoordinates(sum.peak[pol][peak].phi, sum.peak[pol][peak].theta, &RA_CR, &dec_CR); 
    cohSNR = sum.coherent[pol][peak].snr; 
    mcE = sum.mc.energy; 

    getRADec(RA,dec, &sum, &gps, pol, peak); 
    D_grb =doInterpolate(h_grb, RA,dec,unixTime);
    D_grb_12 =doInterpolate(h_grb_12, RA,dec,unixTime);
    D_grb_24 =doInterpolate(h_grb_24, RA,dec,unixTime);
    D_txs = doInterpolate(h_txs,RA,dec,unixTime,true); 
    D_fava = doInterpolate(h_fava,RA,dec,unixTime,true); 
    D_sn = doInterpolate(h_sn,RA,dec,unixTime,true); 
    //now do the offsets

    printf("%d %g %g\n", i, RA,dec); 
    for (int ioffset = 0; ioffset < noffsets; ioffset++) 
    {
      double toffset = gRandom->Uniform(offset_min, offset_max); 
      if (gRandom->Uniform(0,1) < 0.5) toffset *= -1; 
      offsets[ioffset] = toffset; 
      double fakeTime = unixTime + toffset;
      UsefulAdu5Pat fakeGps(d->gps());
      fakeGps.realTime = fakeTime;
      double fakeRA, fakeDec;
      getRADec(fakeRA,fakeDec, &sum, &fakeGps, pol, peak); 
      D_grb_offset[ioffset] =doInterpolate(h_grb, fakeRA,fakeDec,fakeTime);
      D_grb_12_offset[ioffset] =doInterpolate(h_grb_12, fakeRA,fakeDec,fakeTime);
      D_grb_24_offset[ioffset] =doInterpolate(h_grb_24, fakeRA,fakeDec,fakeTime);
      D_txs_offset[ioffset] = doInterpolate(h_txs,fakeRA,fakeDec,fakeTime,true); 
      D_fava_offset[ioffset] = doInterpolate(h_fava,fakeRA,fakeDec,fakeTime,true); 
      D_sn_offset[ioffset] = doInterpolate(h_sn,fakeRA,fakeDec,fakeTime,true); 
    }


    f.cd(); 
    out->Fill(); 
  }
  f.cd(); 
  out->Write(); 
} 


const int source_n_ra_bins = 720;
const int source_n_dec_bins = 360;
const int raw_source_n_ra_bins = 360;
const int raw_source_n_dec_bins = 180;





UCorrelator::KDE2D * fillTMVARADec(double kde_sigma_x, double kde_sigma_y,  TChain * c, double mintime =0, double maxtime = 0, TH2 * h = 0, TH1* dphi_hist = 0, TH1 * dtheta_hist = 0, 
                                   const char * summary_name = "summary", const char * pat_name = "pat", bool onlyOdd = true, bool onlyVPol = true,  
                                   bool mostImpulsive = true, double mc_phi_cut = 4, double mc_theta_cut = 4)
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
    int pol = mostImpulsive? sum->mostImpulsivePolAsInt() : 1; 
    int idx = mostImpulsive ? sum->mostImpulsiveInd() : 0; 
    bool odd = sum->eventNumber & 1; 
    if (onlyOdd &&!odd) continue; 
    if (onlyVPol && pol == 0) continue; 

    if (fabs(sum->peak[pol][idx].theta-sum->mc.theta) < mc_theta_cut && fabs(FFTtools::wrap(sum->peak[pol][idx].phi-sum->mc.phi,360,0)) < mc_phi_cut)
    {
       double ra,dec; 
       UsefulAdu5Pat gps(pat); 
       getRADec(ra,dec, sum, &gps, pol, idx, dtheta_hist, dphi_hist); 
       printf("%d %g %g %g %g %g\n",i, ra,dec, sum->peak[pol][idx].phi, sum->peak[pol][idx].theta, sum->mc.weight); 
       if (h) h->Fill(ra,dec, sum->mc.weight); 
       ras.push_back(ra); 
       decs.push_back(dec); 
       ws.push_back(sum->mc.weight); 
    }
    else
    {
      printf("Skipped event with weight %g\n", sum->mc.weight); 
    }
  }


  delete sum; 
  delete pat; 
  UCorrelator::KDE2D::KDE2DOptions opt(0,kde_sigma_x,kde_sigma_y); 
  return new UCorrelator::KDE2D(ras.size(), &ras[0], &decs[0], &ws[0], opt); 
}




void remakeSourceDistanceHist(const char * key )
{

  TFile f(the_source_distances_file,"UPDATE"); 
  TH3 * h_int = (TH3*) f.Get(Form("%s_int", key)); 
  TH3 * h = (TH3*) UCorrelator::image::makePctileHist(h_int); 
  h->SetName(key); 
  h->SetTitle(key); 
  h->GetXaxis()->SetTitle("RA");
  h->GetYaxis()->SetTitle("dec");
  h->GetZaxis()->SetTitle("time");
  h->Write(key, TObject::kOverwrite); 
}

void makeSourceDistanceHists(const char * dir, const char * key, 
                             const char * time_file, int start_run, int end_run,
                             double kde_sigma_x = 0.5, 
                             double kde_sigma_y = 2, 
                             bool exposure_weighted = true
                             ) 
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
  std::vector<double> raw_ra_bins(raw_source_n_ra_bins+1); 
  std::vector<double> raw_dec_bins(raw_source_n_dec_bins+1); 

  double d_ra = 24./source_n_ra_bins;
  double d_dec = 180./source_n_dec_bins;
  double raw_d_ra = 24./raw_source_n_ra_bins;
  double raw_d_dec = 180./raw_source_n_dec_bins;

  for (int i = 0; i <= source_n_ra_bins; i++) ra_bins[i] = d_ra*i;
  for (int i = 0; i <= source_n_dec_bins; i++) dec_bins[i] = d_dec*i-90;

  for (int i = 0; i <= raw_source_n_ra_bins; i++) raw_ra_bins[i] = raw_d_ra*i;
  for (int i = 0; i <= raw_source_n_dec_bins; i++) raw_dec_bins[i] = raw_d_dec*i-90;


  TFile f(the_source_distances_file,"UPDATE"); 
  TH3 * h = exposure_weighted ? 0 : new TH3F(key,key, source_n_ra_bins, &ra_bins[0], source_n_dec_bins, &dec_bins[0], time_bins.size()-1, &time_bins[0]); 
  TH3 * h_raw = exposure_weighted ? new TH3F(Form("%s_raw",key),key, raw_source_n_ra_bins, &raw_ra_bins[0], raw_source_n_dec_bins, &raw_dec_bins[0], time_bins.size()-1, &time_bins[0]) : 0; 
  TH3 * h_int = exposure_weighted ? new TH3F(Form("%s_int",key),key, source_n_ra_bins, &ra_bins[0], source_n_dec_bins, &dec_bins[0], time_bins.size()-1, &time_bins[0]) : 0; 

  //ok now get the kde's 

  if (h_raw) 
  {
    h_int->SetEntries(h_int->GetNbinsX()*h_int->GetNbinsY()*h_int->GetNbinsZ()); 
    h_raw->SetEntries(h_int->GetNbinsX()*h_int->GetNbinsY()*h_int->GetNbinsZ()); 
  }
  else
  {
    h->SetEntries(h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ()); 
  }

  for (unsigned ti = 1; ti < time_bins.size(); ti++) 
  {
    TH2D this_hist(TString::Format("hradec_%d",ti), TString::Format("RA dec %d", ti), raw_source_n_ra_bins, 0 , 24, raw_source_n_dec_bins,-90,90); 
    UCorrelator::KDE2D * kde = fillTMVARADec(kde_sigma_x, kde_sigma_y,  &c, time_bins[ti-1], time_bins[ti], &this_hist); 

    if (kde->getN() == 0) 
    {
      printf("No events in time slice!"); 

      if (!exposure_weighted) 
      {
        for (int yi = 1; yi <= source_n_dec_bins; yi++) 
        {
          for (int xi = 1; xi<= source_n_ra_bins; xi++) 
          {
            h->SetBinContent(xi,yi,ti, 1);
          }
        }
      }
      delete kde;
      continue;
    }
    else if (exposure_weighted) 
    {
      for (int yi = 1; yi <= h_raw->GetNbinsY(); yi++) 
      {
        for (int xi = 1; xi<= h_raw->GetNbinsX(); xi++) 
        {
          h_raw->SetBinContent(xi,yi,ti, this_hist.GetBinContent(xi,yi)); 
        }
      }
    }

    //now let's fill the  histogram! 

    TH2D htmp("htmp", "htmp", source_n_ra_bins, &ra_bins[0], source_n_dec_bins, &dec_bins[0]);

    for (int yi = 1; yi <= htmp.GetNbinsY(); yi++) 
    {
      for (int xi = 1; xi <= htmp.GetNbinsX(); xi++) 
      {
        double ra = htmp.GetXaxis()->GetBinCenter(xi); 
        double dec = htmp.GetYaxis()->GetBinCenter(yi); 
       //get wrap-around right! 
        double val = kde->operator()(ra,dec) + kde->operator()(ra+24,dec) + kde->operator()(ra-24,dec);
        htmp.SetBinContent(xi,yi, val);
        if (exposure_weighted) 
        {
          h_int->SetBinContent(xi,yi,ti,val); 
        }
      }
    }

    delete kde; 

    if (!exposure_weighted) 
    {
      TCanvas * can = new TCanvas; 
      can->Divide(3,1);
      can->cd(1);
      htmp.DrawCopy("col2z"); 
      TH2 * pct = (TH2*) UCorrelator::image::makePctileHist(&htmp); 
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
  }

  f.cd();
  if (exposure_weighted) 
  {
    h = (TH3*) UCorrelator::image::makePctileHist(h_int); 
  }

  h->SetName(key); 
  h->SetTitle(key); 

  h->GetXaxis()->SetTitle("RA");
  h->GetYaxis()->SetTitle("dec");
  h->GetZaxis()->SetTitle("time");

  h_raw->GetXaxis()->SetTitle("RA");
  h_raw->GetYaxis()->SetTitle("dec");
  h_raw->GetZaxis()->SetTitle("time");

  h_int->GetXaxis()->SetTitle("RA");
  h_int->GetYaxis()->SetTitle("dec");
  h_int->GetZaxis()->SetTitle("time");


  h->Write(key, TObject::kOverwrite); 
  if (exposure_weighted) 
  {
    h_raw->Write(h_raw->GetName(), TObject::kOverwrite); 
    h_int->Write(h_int->GetName(), TObject::kOverwrite); 
  }
}

void allSourceDistanceHists(int which = 15) 
{

  if (which & 1) 
  {
    makeSourceDistanceHists("simulated_txs_src","txs",0,1,100, 0.125,1);
  }
  if (which & 2) 
  {
    makeSourceDistanceHists("simulated_grb_src","grb","/project2/kicp/cozzyd/mc/icemc_grb/src/run1/source_times.txt",1,300);
    makeSourceDistanceHists("simulated_grb_src_12","grb_12","/project2/kicp/cozzyd/mc/icemc_grb/src_12/run1/source_times.txt",1,100);
    makeSourceDistanceHists("simulated_grb_src_24","grb_24","/project2/kicp/cozzyd/mc/icemc_grb/src_24/run1/source_times.txt",1,100);
  }
  if (which & 4) 
  {
    makeSourceDistanceHists("simulated_sn_src","sn","/project2/kicp/cozzyd/mc/icemc_sn/src/run1/source_times.txt",1,100);
  }
  if (which & 8) 
  {
    makeSourceDistanceHists("simulated_fava_src","fava","/project2/kicp/cozzyd/mc/icemc_fava/src/run1/source_times.txt",1,200);
  }

}





void allOverlaps(int which =  31)
{
  if (which & 1) 
  {
    fillOverlapHists("source_maps_eval/*full.root","A3Overlaps.root",0) ;
  }

  if (which & 2) 
  {
    setenv("ANITA_MC_DATA","/project2/kicp/cozzyd/mc/icemc_fava/src",1);
    fillOverlapHists("source_maps_eval/*simulated_fava_src.root","fava_overlaps.root","simulated_fava_src");
  }

  if (which & 4) 
  {
    setenv("ANITA_MC_DATA","/project2/kicp/cozzyd/mc/icemc_sn/src",1);
    fillOverlapHists("source_maps_eval/*simulated_sn_src.root","sn_overlaps.root","simulated_sn_src");
  }



  if (which & 8) 
  {
    setenv("ANITA_MC_DATA","/project2/kicp/cozzyd/mc/icemc_txs/src",1);
    fillOverlapHists("source_maps_eval/*simulated_txs_src.root","txs_overlaps.root","simulated_txs_src");
  }


  if (which & 16) 
  {
    setenv("ANITA_MC_DATA","/project2/kicp/cozzyd/mc/icemc_grb/src",1);
    fillOverlapHists("source_maps_eval/*simulated_grb_src.root","grb_overlaps.root","simulated_grb_src");

    setenv("ANITA_MC_DATA","/project2/kicp/cozzyd/mc/icemc_grb/src_12",1);
    fillOverlapHists("source_maps_eval/*simulated_grb_src_12.root","grb_12_overlaps.root","simulated_grb_src_12");

    setenv("ANITA_MC_DATA","/project2/kicp/cozzyd/mc/icemc_grb/src_24",1);
    fillOverlapHists("source_maps_eval/*simulated_grb_src_24.root","grb_24_overlaps.root","simulated_grb_src_24");
  }
}



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

int nbinsD = 100; 
int nbinsF = 30; 
double minF = 0;
double maxF = 3; 

int nbinslogO =70; 
double minlogO = -1; 
double maxlogO = 6; 


TH3D * getEfficiency(TTree * overlap_tree, const char * thermal_tree_file,  const char * key, double fraction = 1, bool onlyEven = true) 
{

   
  //get the denominator from the thermal tree
 

  TFile fth(thermal_tree_file); 
  TTree * th = (TTree*) fth.Get("simulation"); 

  th->Draw("F", onlyEven ? "weight * isMostImpulsive * ((eventNumber_nnn & 1) == 0)": "weight * isMostImpulsive","goff"); 

  TH1 * htemp = (TH1*) gDirectory->Get("htemp");
  double denom = htemp->Integral()*fraction; 
  printf("Denom: %g\n", denom); 

  fth.Close(); 
  gROOT->cd(); 

  int N = overlap_tree->Draw(TString::Format("F:O:D_%s",key), onlyEven ? "weight * ( (event & 1)==0) * (pol == 1) ": "weight","goff"); 

  double *W = overlap_tree->GetW(); 
  double *vF = overlap_tree->GetV1();
  double *vO = overlap_tree->GetV2();
  double *vD = overlap_tree->GetV3();

  TH3D * eff = new TH3D(TString::Format("%s__eff", key), TString::Format("Efficiency Estimate for %s", key), 
                       nbinsF,minF, maxF,  //F
                       nbinslogO, minlogO, maxlogO,  //log O 
                       nbinsD, 0, 1);  // D




  for (int k =1;  k <= eff->GetNbinsZ(); k++) 
  {
    double D =eff->GetZaxis()->GetBinLowEdge(k);
    for (int j = 1; j <= eff->GetNbinsY(); j++) 
    {
      double O = pow(10,-eff->GetYaxis()->GetBinLowEdge(j)); 
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
       // printf("%g %g %g %g\n", D, -log10(O), F, num); 

        eff->SetBinContent(i,j,k,num/denom); 
      }
    }
  }

  //delete tmp; 
  return eff; 

}

double sigmaEff = 0.1; 
std::map<double, UCorrelator::CachedFC*> grids; 

TH3D * getScaledSup(TH3D * eff, TH3D * bg, double CL=0.9, double poisson_prior=1, int N = 100000)
{
  TString keys(bg->GetName()); 
  keys.ReplaceAll("__bg",""); 
  const char * key = keys.Data(); 
  TH3D * scaled = new TH3D(Form("%s_scaled_sup_%g_%g", key, CL*1000, poisson_prior*10), Form("%s,  #bar{s_{up}}/#epsilon, CL=%g, poisson_prior =%g", key, CL, poisson_prior),  
                           bg->GetNbinsX(), bg->GetXaxis()->GetXmin(), bg->GetXaxis()->GetXmax(), 
                           bg->GetNbinsY(), bg->GetYaxis()->GetXmin(), bg->GetYaxis()->GetXmax(), 
                           bg->GetNbinsZ(), bg->GetZaxis()->GetXmin(), bg->GetZaxis()->GetXmax()) ;


  ROOT::Math::GSLRandomEngine gslran; 
  gslran.Initialize(); 

  UCorrelator::CachedFC * grid; 
  if (!grids.count(CL))
  {
    grid = new UCorrelator::CachedFC(CL,350,300,0.25); 
    grids[CL] = grid; 
  }
  else
  {
    grid = grids[CL]; 
  }

  for (int k = 1; k <= bg->GetNbinsZ(); k++) 
  {
    for (int j = 1; j <= bg->GetNbinsY(); j++) 
    {
      for (int i = 1; i <=bg->GetNbinsX(); i++) 
      {
        double bg_point = bg->GetBinContent(i,j,k); 
        double fc_sum = 0; 
        double eff_pt = eff->GetBinContent(i,j,k); 
        int this_N = N; 
        if (bg_point > 8) //probably not using this, so don't try so hard.  
        {
          this_N /= (bg_point*bg_point); 
        }
        if (eff_pt < 0.5) //probably not using this, so don't try so hard
        {
          this_N *= 4*eff_pt * eff_pt; 
        }

        if (this_N < 10) this_N = 10; 

        for (int t = 0; t < this_N; t++) 
        {
          double eff = eff_pt; 
          if (sigmaEff) 
          {
            do
            {
              eff = gRandom->Gaus(eff_pt, sigmaEff); 
            } while ( eff > 1 || eff < 0); 
          }

          //interpretation: we have a point estimate of the background in a sideband that
          //is 7 times bigger. 
          double bg_realization = gslran.Gamma(poisson_prior+7*bg_point, 1./7); 
          double seen = gRandom->Poisson(bg_realization); 
          if (seen > 300 || bg_realization > 250) 
          {
            printf("OOB: seen: %f, bg: %f\n", seen, bg_realization); 

          }
          double fc_val = grid->upperLimit(seen,bg_realization); 
          fc_sum += fc_val/eff; 
        }

        double fc_mean = fc_sum/this_N; 
        if (fc_mean == 0) fc_mean = 100; //hack! 
        scaled->SetBinContent(i,j,k, fc_mean); 
        printf("%d %d %d %g %g %g\n", i,j,k, bg_point, eff_pt, fc_mean);
      }
    }
  }

  return scaled; 
}



TH3D * getBackgroundEstimate(TTree * overlap_tree, const char * key) 
{
  TH3D * bg = new TH3D(TString::Format("%s__bg", key), TString::Format("BG Estimate for %s;F;-log(O);D", key), 
                       nbinsF, minF, maxF,  //F
                       nbinslogO, minlogO, maxlogO,  //log O 
                       nbinsD, 0, 1);  // D


  for (int k =1;  k <= bg->GetNbinsZ(); k++) 
  {
    double D = bg->GetZaxis()->GetBinLowEdge(k); 
    TH2 * tmp = new TH2D(TString::Format("%s_tmp_bg", key), TString::Format("Differential BG Estimate for %s", key), 
                       nbinsF, minF, maxF,  //F
                       nbinslogO, minlogO, 12);  //log O 

    overlap_tree->Draw(Form("-log10(O+1e-15):F >> %s_tmp_bg", key), Form("(pol == 1) * (D_%s_offset < %g)", key, D),"goff");
    for (int j = 1; j <= bg->GetNbinsY(); j++) 
    {
      double nlogO = bg->GetYaxis()->GetBinLowEdge(j); 
      for (int i = 1; i <= bg->GetNbinsX(); i++) 
      {
        double F = bg->GetXaxis()->GetBinLowEdge(i); 
        double sum = getUpperRightSum(tmp,i,j,true); 
        double scale = 0.01; 
//        printf(" D=%g -log O =%g  F=%g sum = %g  scale =%g\n",D, nlogO,F,  sum, scale); 
        bg->SetBinContent(i,j,k, sum*scale); //jefrey's prior 
        bg->SetBinError(i,j,k, (0.5+sqrt(sum+0.25))*scale);
      }
    }
    delete tmp; 
  }

  return bg; 

}

/*
TH3D * getOldBackgroundEstimate(TTree * overlap_tree, TH3 * source_distances, const char * key) 
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
*/


void optimization(const char * overlap_file, const char * out_file, const char * key, const char * thermal_trees)
{
  TFile f_sim_overlaps(overlap_file); 
  TTree * overlap = (TTree*)f_sim_overlaps.Get("overlap"); 

  TFile f_data_overlaps("A3Overlaps.root"); 
  TTree * data_overlap = (TTree*) f_data_overlaps.Get("overlap"); 

  TFile fout(out_file,"RECREATE"); 
  TH3D * eff = getEfficiency(overlap,thermal_trees,key); 
  fout.cd(); 
  eff->Write();
  TH3D * bg = getBackgroundEstimate(data_overlap,key); 
  fout.cd(); 
  bg->Write(); 

//  TH3* scaled_sup_90 = getScaledSup(eff,bg,0.9,0.5); 
//  TH3* scaled_sup_975 = getScaledSup(eff,bg,0.975,0.5); 

  TH3* scaled_sup_90_uniform = getScaledSup(eff,bg,0.9,1); 
//  TH3* scaled_sup_975_uniform = getScaledSup(eff,bg,0.975,1); 
  fout.cd(); 
 // scaled_sup_90->Write();
 // scaled_sup_975->Write();
  scaled_sup_90_uniform->Write();
 // scaled_sup_975_uniform->Write();

}

void favaOptimization()
{
  optimization("fava_overlaps.root","fava_opt.root","fava","thermalTrees/simulated_fava_src_1-200_sinsub_10_3_ad_2.root");
}

void snOptimization()
{
  optimization("sn_overlaps.root","sn_opt.root","sn","thermalTrees/simulated_sn_src_1-200_sinsub_10_3_ad_2.root");
}

void grbOptimization()
{
  optimization("grb_overlaps.root","grb_opt.root","grb","thermalTrees/simulated_grb_src_1-300_sinsub_10_3_ad_2.root");
}


void grb12Optimization()
{
  optimization("grb_12_overlaps.root","grb_12_opt.root","grb_12","thermalTrees/simulated_grb_src_12_1-100_sinsub_10_3_ad_2.root");
}

void grb24Optimization()
{
  optimization("grb_24_overlaps.root","grb_24_opt.root","grb_24","thermalTrees/simulated_grb_src_24_1-100_sinsub_10_3_ad_2.root");
}




void txsOptimization()
{


  optimization("txs_overlaps.root","txs_opt.root","txs","thermalTrees/simulated_txs_src_1-100_sinsub_10_3_ad_2.root");

  /*
  TFile f_txs_overlaps("txs_overlaps.root"); 
  TTree * overlap = (TTree*)f_txs_overlaps.Get("overlap"); 

  TFile f_data_overlaps("A3Overlaps.root"); 
  TTree * data_overlap = (TTree*) f_data_overlaps.Get("overlap"); 

  TFile fout("txs_opt.root","RECREATE"); 
  TH3D * eff = getEfficiency(overlap,"thermalTrees/simulated_txs_src_1-100_sinsub_10_3_ad_2.root","txs"); 
  fout.cd(); 
  eff->Write();
  TH3D * bg = getBackgroundEstimate(data_overlap,"txs"); 
  fout.cd(); 
  bg->Write(); 

  TH3* scaled_sup_90 = getScaledSup(eff,bg,0.9); 
  TH3* scaled_sup_975 = getScaledSup(eff,bg,0.975); 

  fout.cd(); 
  scaled_sup_90->Write();
  scaled_sup_975->Write();
  */

}

void allOptimization() 
{

  favaOptimization(); 
  txsOptimization(); 
  snOptimization(); 
  grbOptimization(); 
  grb12Optimization();
  grb24Optimization(); 
}
