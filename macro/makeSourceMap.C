#include "FFTtools.h" 
#include "AnitaConventions.h" 
#include "AnitaDataset.h" 

bool isMC = 1;

double cutoff = 3; 
const char * weight = "F > 3";
// const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && theta < 40 )";
// const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && isMostImpulsive && !payloadBlast && MaxPeak < 1000 && theta < 40 && ( (HPolTrigger && iteration < 5) || (VPolTrigger && iteration > 4))  && !isPulser  )";

UCorrelator::ProbabilityMap::Params * map_params()
{

  StereographicGrid * g= new StereographicGrid(2024,2024); 

  TF1 * f_dtheta = new TF1("ftheta", "[0] / x^[1]", 1, 50);
  TF1 * f_dphi = new TF1("fphi", "[0] / x^[1]", 1, 50);
  //anita4 fit from wais.
  f_dtheta->SetParameter(0, 5.431); 
  f_dtheta->SetParameter(1, 1.155); 
  f_dphi->SetParameter(0, 4.061); 
  f_dphi->SetParameter(1, 1.038);
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
  p->seg = g; 
  p->point = resolutionModel; 
  p->collision_detection = false; 
 

  return p; 

}



void addRuns(TChain & c, int start_run, int end_run)
{
  for (int i = start_run; i <= end_run; i+=40) 
  {
    if (start_run < i + 40 && end_run >= i)
    {
      TString adding;
      if(isMC == 0){
        adding = TString::Format("thermalTrees/a4all_%d-%d_max_30001_sinsub_10_3_ad_2.root",i,i+39);
      }else{
        adding = TString::Format("thermalTrees/simulated_%d-%d_max_1000_sinsub_10_3_ad_2.root",i,i+39);
      }
      printf("Adding %s\n", adding.Data() ); 
      c.Add(adding.Data() ); 
    }
  }
}

std::set<int> * getRemovedEvents(const char * file, std::vector<int>  * runs = 0, std::vector<int> * events = 0, std::vector<int> * iters = 0) 
{
  FILE * f = fopen(file,"r") ; 
  if (!f) return 0; 
  std::set<int> * removed = new std::set<int>; 
  char buf[1024]; 
  while (fgets(buf,sizeof(buf),f))
  {
    char * comment =strchr(buf,'#'); 
    if (comment) *comment=0; 

    int run,event,iteration; 
    sscanf(buf,"%d %d %d", &run,&event,&iteration); 

    if (runs) runs->push_back(run); 
    if (events) events->push_back(event); 
    if (iters) iters->push_back(iteration); 
    removed->insert(event); 
  }

  return removed; 
}

// int makeSourceMap(int start_run = 50, int end_run = 367, const char * prefix = "source_maps/test_")
int makeSourceMap(int start_run = isMC?1:50, int end_run = isMC?200:367, const char * prefix = "source_maps/test_")
{

  // Start getting the run / event numbers of events that pass our cut

  TChain c(isMC?"simulation":"anita4"); 

  addRuns(c,start_run,end_run); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  int total_event_n = c.Draw("run:eventNumber:iteration:F",  TCut(TString::Format("(%s) && (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 
  printf("%d events pass selection\n", total_event_n); 


  UCorrelator::ProbabilityMap::Params * p = map_params(); 

  TFile * f  = new TFile(TString::Format("%s%d_%d.root",prefix,start_run, end_run), "RECREATE"); 
  // UCorrelator::ProbabilityMap *map_weighted = new UCorrelator::ProbabilityMap(p); 
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(p); 

  TTree * tr = new TTree("events","events"); 
  int run, ev, pol, peak,  nsegs ; 
  double S, F, dinteg, dinteg_norm; 

  // tr->Branch("event",&ev); 
  // tr->Branch("run",&run); 
  // tr->Branch("pol",&pol); 
  // tr->Branch("peak",&peak); 
  // tr->Branch("S",&S); 
  // tr->Branch("F",&F); 
  // tr->Branch("dinteg",&dinteg); 
  // tr->Branch("dinteg_norm",&dinteg_norm); 
  // tr->Branch("nsegs",&nsegs); 
  

  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  double last_integ = 0; 
  double last_integ_norm = 0; 
  for (int i = 0; i < total_event_n; i+=1) 
  {
    run = c.GetV1()[i]; 
    if (run!= loaded_run)
    {
      if (sumfile) delete sumfile;
      if(isMC == 0){
        sumfile = new TFile(TString::Format("/Volumes/SDCard/data/%s/%d_max_30001_sinsub_10_3_ad_2.root", "a4all",run));
      }else{
        sumfile = new TFile(TString::Format("/Volumes/SDCard/data/%s/%d_max_1000_sinsub_10_3_ad_2.root", "simulated",run));
      } 
         
      gROOT->cd(); 
      sumtree = (TTree*) sumfile->Get(isMC?"simulation":"anita4"); 
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
    double integ = map->getProbSumsIntegral(false); 
    double integ_norm = map->getProbSumsIntegral(true); 
    dinteg = integ-last_integ; 
    dinteg_norm = integ_norm-last_integ_norm; 
    last_integ_norm = integ_norm; 
    last_integ = integ_norm; 

    // tr->Fill(); 

    if (F > cutoff) 
    {
      map->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, 1); 
    }


  }


  f->cd(); 
  map->Write("map_unweighted"); 
  // map_weighted->Write("map_weighted"); 
  tr->Write(); 
  return 0; 
}




