#include "FFTtools.h" 
#include "AnitaConventions.h" 
#include "AnitaDataset.h" 

double cutoff = 1; 
const char * weight = "F > 1";
// const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && theta < 40 )";
// const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && isMostImpulsive && !payloadBlast && MaxPeak < 1000 && theta < 40 && ( (HPolTrigger && iteration < 5) || (VPolTrigger && iteration > 4))  && !isPulser  )";

UCorrelator::ProbabilityMap::Params * map_params()
{
  // pixel x, pixel y , max meter x, max meter y
  StereographicGrid * g= new StereographicGrid(1000,1000,2000000,2000000); 

  TF1 * f_dtheta = new TF1("ftheta", "[0] / x^[1]", 1, 50);
  TF1 * f_dphi = new TF1("fphi", "[0] / x^[1]", 1, 50);
  //anita4 fit from wais.
  f_dtheta->SetParameter(0, 5.431); 
  f_dtheta->SetParameter(1, 1.155); 
  f_dphi->SetParameter(0, 28.87); 
  f_dphi->SetParameter(1, 1.398);
  //anita3
  // f_dtheta->SetParameter(0, 0.3936); 
  // f_dtheta->SetParameter(1, 0.2102); 
  // f_dphi->SetParameter(0, 1.065); 
  // f_dphi->SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel * snrResolutionModel = new UCorrelator::PointingResolutionParSNRModel (*f_dtheta, *f_dphi, true,true);
  UCorrelator::PointingResolutionModelPlusHeadingError * resolutionModel = new UCorrelator::PointingResolutionModelPlusHeadingError(20, snrResolutionModel); 

  Refraction::SphRay * ref = new Refraction::SphRay; 

  UCorrelator::ProbabilityMap::Params *p = new UCorrelator::ProbabilityMap::Params; 
  p->refract = ref; 
  // p->refract = 0; 
  p->seg = g; 
  p->point = resolutionModel; 
  p->collision_detection = true; 
  p->verbosity = 0; // verbosity level for output info.
  p->maximum_distance = 2.5;
 

  return p; 

}



void addRuns(TChain & c, int start_run, int end_run, const char* thermalTreeFormat)
{
  for (int i = start_run; i <= end_run; i+=40) 
  {
    if (start_run < i + 40 && end_run >= i)
    {
      TString adding = TString::Format(thermalTreeFormat,i,i+39);
      printf("Adding %s\n", adding.Data() ); 
      c.Add(adding.Data() ); 
    }
  }
}

int _makeSourceMap(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run = 50, int end_run =367, const char * outputFilePrefix = "source_maps/test_")
{

  // Start getting the run / event numbers of events that pass our cut

  TChain c(treeName); 

  addRuns(c,start_run,end_run, thermalTreeFormat); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  int total_event_n = c.Draw("run:eventNumber:pol*5+peak:F",  TCut(TString::Format("(%s) && (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 
  printf("%d events pass selection\n", total_event_n); 


  UCorrelator::ProbabilityMap::Params * p = map_params(); 

  TFile * f  = new TFile(TString::Format("%s%d_%d.root",outputFilePrefix,start_run, end_run), "RECREATE"); 
  // UCorrelator::ProbabilityMap *map_weighted = new UCorrelator::ProbabilityMap(p); 
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(p); 

  TTree * tr = new TTree("events","events"); 
  int run, ev, pol, peak,  nsegs ; 
  double S, F, dinteg, dinteg_norm; 

  tr->Branch("event",&ev); 
  tr->Branch("run",&run); 
  tr->Branch("pol",&pol); 
  tr->Branch("peak",&peak); 
  tr->Branch("S",&S); 
  tr->Branch("F",&F); 
  tr->Branch("dinteg",&dinteg); 
  tr->Branch("dinteg_norm",&dinteg_norm); 
  tr->Branch("nsegs",&nsegs);
  

  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  double last_integ = 0; 
  double last_integ_norm = 0; 
  // for (int i = 0; i < total_event_n; i+=1) 
  for (int i = 0; i < total_event_n; i+=1) 
  {
    run = c.GetV1()[i];
    // skip 45 degree wais runs. 
    // if (run<136){
    //   continue;
    // }
    if (run!= loaded_run)
    {
      if (sumfile) delete sumfile;
      sumfile = new TFile(TString::Format(summaryFileFormat,run));
         
      gROOT->cd(); 
      sumtree = (TTree*) sumfile->Get(treeName); 
      sumtree->SetBranchAddress("summary",&sum); 
      sumtree->SetBranchAddress("pat",&gps); 
      sumtree->BuildIndex("eventNumber"); 
      loaded_run = run; 
    }
      
    ev = int(c.GetV2()[i]); 
    S = c.GetW()[i]; 
    sumtree->GetEntryWithIndex(ev); 
    pol = int(c.GetV3()[i]) / 5; 
    peak = int(c.GetV3()[i]) % 5; 
    F = c.GetV4()[i]; 
    printf("index = %d \t run = %d \t eventNumber = %d \t F = %g \t S = %g \n",i,run,ev,F,S); 

    // nsegs = map_weighted->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, S); 
    nsegs = map->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, S);
    std::cout<< "\tsnr = "<< sum->deconvolved_filtered[pol][peak].snr << " longitude="<<sum->peak[pol][peak].longitude<<" latitude"<<sum->peak[pol][peak].latitude<< std::endl; 
    double integ = map->getProbSumsIntegral(false); 
    double integ_norm = map->getProbSumsIntegral(true); 
    dinteg = integ-last_integ; 
    dinteg_norm = integ_norm-last_integ_norm; 
    last_integ_norm = integ_norm; 
    last_integ = integ_norm; 

    tr->Fill(); 

    // if (F > cutoff) 
    // {
    //   map->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, 1); 
    // }
  }
  f->cd(); 
  map->Write("map_unweighted"); 
  // map_weighted->Write("map_weighted"); 
  tr->Write(); 
  return 0; 
}

void makeSourceMap(const char * treeName){
  int start_run,end_run;
  if (!strcmp(treeName,"wais")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/wais/%d_max_30001_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/wais_%d-%d_max_30001_sinsub_10_3_ad_2.root";
    start_run = 120;
    end_run = 155;
    _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run);

  }else if(!strcmp(treeName,"anita4")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/a4all/%d_max_30001_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/a4all_%d-%d_max_30001_sinsub_10_3_ad_2.root";
    start_run = 50;
    end_run = 367;
    _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run);
  }else if(!strcmp(treeName,"simulation")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/simulated/%d_max_501_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/simulated_%d-%d_max_501_sinsub_10_3_ad_2.root";
    start_run = 1;
    end_run = 10;
    _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run);
  }else{
    std::cout<< "wrong input treeName"<<std::endl;
  }
}




