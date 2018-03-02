#include "FFTtools.h" 
#include "AnitaConventions.h" 
#include "AnitaDataset.h" 



double cutoff = 1; 
const char * weight = "F>1";





UCorrelator::ProbabilityMap::Params * map_params()
{

  StereographicGrid * g= new StereographicGrid(1024,1024); 

  TF1 * f_dtheta = new TF1("ftheta", "[0] / x^[1]", 1, 50);
  TF1 * f_dphi = new TF1("fphi", "[0] / x^[1]", 1, 50);
  f_dtheta->SetParameter(0, 0.3936); 
  f_dtheta->SetParameter(1, 0.2102); 
  f_dphi->SetParameter(0, 1.065); 
  f_dphi->SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel * m1 = new UCorrelator::PointingResolutionParSNRModel (*f_dtheta, *f_dphi, true,true);
  UCorrelator::PointingResolutionModelPlusHeadingError * m = new UCorrelator::PointingResolutionModelPlusHeadingError(20, m1); 

  Refraction::SphRay * ref = new Refraction::SphRay; 

  UCorrelator::ProbabilityMap::Params *p = new UCorrelator::ProbabilityMap::Params; 
  p->refract = ref; 
  p->seg = g; 
  p->point = m; 
  p->collision_detection = false; 
 

  return p; 

}


//add the thermalTrees
void addRuns(TChain & c, int start_run, int end_run)
{
  for (int i = start_run; i <= end_run; i+=40) 
  {
    if (start_run < i + 40 && end_run >= i)
    {
      TString adding = TString::Format("thermalTrees/a4all_%d-%d_max_30001_sinsub_10_3_ad_2.root",i,i+39);
      // TString adding = TString::Format("thermalTrees/simulated_%d-%d_max_1000_sinsub_10_3_ad_2.root",i,i+39);
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

    int run,eventNumber,iteration; 
    sscanf(buf,"%d %d %d", &run,&eventNumber,&iteration); 

    if (runs) runs->push_back(run); 
    if (events) events->push_back(eventNumber); 
    if (iters) iters->push_back(iteration); 
    removed->insert(eventNumber); 
  }

  return removed; 
}

//input map_weighted and output map_weighted_culled, removed some events from list.
int removeEvents(const char * maps_file = "all_source_maps.root",
                 const char * removed_file = "removed_events.txt",
                 const char * in_key = "map_weighted", const char * out_key =  "map_weighted_culled") 
{
  std::vector<int> runs; 
  std::vector<int> events; 
  std::vector<int> iters; 

  getRemovedEvents(removed_file, &runs,&events,&iters); 

  TFile f(maps_file,"UPDATE"); 
  UCorrelator::ProbabilityMap * m = (UCorrelator::ProbabilityMap*) f.Get(in_key); 
  printf("%x\n", m); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  TTree * eventstree = 0; 
  TFile * eventsfile = 0; 

  double S; 
  int loaded_run = 0; 
  for (unsigned i = 0; i < events.size(); i++) 
  {
    int run = runs[i]; 

    if (run!=loaded_run) 
    {
      if (sumfile) delete sumfile; 
      sumfile = new TFile(TString::Format("/Volumes/SDCard/data/a4all/%d_max_30001_sinsub_10_3_ad_2.root",run));; 
      gROOT->cd(); 

      if (eventsfile) delete eventsfile; 
      eventsfile = new TFile(TString::Format("source_maps/%d_%d.root", run,run)); 
      gROOT->cd(); 

      sumtree = (TTree*) sumfile->Get("anita4"); 
      eventstree = (TTree*) eventsfile->Get("events"); 

      sumtree->SetBranchAddress("summary",&sum); 
      sumtree->SetBranchAddress("pat",&gps); 
      sumtree->BuildIndex("eventNumber"); 
      eventstree->BuildIndex("eventNumber"); 
      eventstree->SetBranchAddress("S",&S); 

      loaded_run =run; 
    }

    int eventNumber = events[i]; 
    int iteration = iters[i]; 
    sumtree->GetEntryWithIndex(eventNumber); 
    eventstree->GetEntryWithIndex(eventNumber); 
    printf("%d %d %d %g\n", run,eventNumber, iteration, S); 
    m->add(sum,gps,AnitaPol::AnitaPol_t (iteration / 5), iteration % 5, -S); 
  }
  
  f.cd(); 
  m->Write(out_key); 
  return events.size(); 
}


int evaluateSourceMap(int start_run = 50, int end_run = 367, 
                      bool mc = false, const char * prefix = "source_maps_eval/",
                      const char * sourceMapName="source_maps/50_367.root",
                      const char * sourceMapTree = "map_weighted", 
                      // const char * removed_events_file = "removed_events.txt", 
                      const char * removed_events_file = NULL, 
                      bool weighted = true, bool vpol_only = false ) 
{
  //Chain through the thermal tree.(tmva feature tree with Fisher value)
  TChain thermalChain(mc ? "simulation":"anita4"); 
  if (mc) {
//    thermalChain.Add("thermalTrees/simulated*.root"); 
    thermalChain.Add("thermalTrees/simulated*.root"); 
  }else{
    addRuns(thermalChain,start_run,end_run); 
  }

  TFile outputFile(TString::Format("%s%d_%d_%s.root",prefix,start_run,end_run, mc ? "sim" : "full" ),"RECREATE"); 
  //output file and output tree named overlap
  TTree * outputTree = new TTree("overlap","Overlap"); 
  double O; 
  double S; 
  int run;
  int eventNumber;
  int pol; 
  double theta; 
  double base_sum; 
  double F; 
  double polangle; 
  double wgt = 1; 
  double mcE = 0; 
  int peak = 0;
  int nclustered[10]; 
  int max_base_index; 
  double max_base_p; 
  int removed = 0; 

  //branches in outputTree, the main purpose are to fill those branches.
  outputTree->Branch("O",&O); 
  outputTree->Branch("S",&S); 
  outputTree->Branch("run",&run); 
  outputTree->Branch("F",&F); 
  outputTree->Branch("eventNumber",&eventNumber); 
  outputTree->Branch("pol",&pol); 
  outputTree->Branch("theta",&theta); 
  outputTree->Branch("nclustered",&nclustered,"nclustered[10]/I"); 
  outputTree->Branch("polangle",&polangle); 
  outputTree->Branch("base_sum",&base_sum); //ps sum from all bases.
  outputTree->Branch("max_base_index",&max_base_index); // max ps base's id
  outputTree->Branch("max_base_p",&max_base_p); //max ps base's ps
  outputTree->Branch("weight",&wgt); // simulation event's weight 
  outputTree->Branch("mcE",&mcE); // simlation evetns's energy
  outputTree->Branch("peak",&peak); 
  outputTree->Branch("removed",&removed); // is this event removed or not.

  TFile sourceMapFile(sourceMapName); 
  UCorrelator::ProbabilityMap * map = (UCorrelator::ProbabilityMap*) sourceMapFile.Get(sourceMapTree); 
  UCorrelator::ProbabilityMap * source_map = map; //for mc 
  UCorrelator::ProbabilityMap::Params * map_pars = map_params(); 
  int total_events_n = thermalChain.Draw("run:eventNumber:iteration:F",  TCut(TString::Format("(%s) * (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 
  std::cout<< "total number of events: "<< total_events_n<< std::endl;
  std::vector<std::vector<double> > counts(10, std::vector<double> (map->segmentationScheme()->NSegments())); 

  for (int d = 0; d < 10; d++){
    // for ten levels, group the adjacent segment wgt. Store into counts.
    map->groupAdjacent(map->getWgtAboveLevel(d), 0, &counts[d][0]); 
  }


  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  //if removed_events_file exsit. otherwise removed_events be 0
  std::set<int> * removed_events = removed_events_file ? getRemovedEvents(removed_events_file) : 0; 

  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  for (int i = 0; i < total_events_n; i++){
    //loop through thermal chain
    run = thermalChain.GetV1()[i]; 
    if (run!= loaded_run)
    {
        if (sumfile) delete sumfile; 
        sumfile = new TFile(TString::Format("/Volumes/SDCard/data/%s/%d_max_30001_sinsub_10_3_ad_2.root", mc ? "simulated" : "a4all",run));; 
        gROOT->cd(); 
        sumtree = (TTree*) sumfile->Get(mc ? "simulation" : "anita4"); 
        sumtree->SetBranchAddress("summary",&sum); 
        sumtree->SetBranchAddress("pat",&gps); 
        sumtree->BuildIndex("eventNumber"); 
        loaded_run =run; 
    }
    //each eventNumber in thermalChain, get S, pol, peak, F, eventNumber
    eventNumber = int(thermalChain.GetV2()[i]); 
    S = thermalChain.GetW()[i]; 
    pol = int(thermalChain.GetV3()[i]) / 5; 
    peak = int(thermalChain.GetV3()[i]) % 5;
    F = thermalChain.GetV4()[i]; 
    //get more variables from summary file
    sumtree->GetEntryWithIndex(eventNumber);  
    polangle = 90/TMath::Pi() * TMath::ATan2(sum->deconvolved_filtered[pol][peak].max_dU, sum->deconvolved_filtered[pol][peak].max_dQ); 
    theta = sum->peak[pol][peak].theta; 

    //this eventNumber is in the removed list.
    removed = (removed_events && removed_events->count(eventNumber));  
    if (mc){
      wgt = sum->mc.weight; 
      mcE = sum->mc.energy; 
      if ( fabs(sum->mc.theta-theta) > 4) continue; 
      if ( fabs(FFTtools::wrap(sum->mc.phi-sum->peak[pol][peak].phi, 360,0)) > 4) continue; 
    }


    if ( (!weighted && F < cutoff) || (vpol_only && pol == 0)) continue; 

    std::vector<std::pair<int,double> > base_contribution;
    //The maximun density of each bin.
    std::vector<std::pair<int,double> > max_dens;
    double inv_two_pi_sqrt_det; 

    if (mc){
      // if mc, create new map, add this current event, combine with the initial source map.
      map = new UCorrelator::ProbabilityMap(map_pars); 
      map->add(sum,gps,AnitaPol::AnitaPol_t(pol),peak,S); 
      map->combineWith(*source_map); 
      for (int d = 0; d< 10; d++){
        map->groupAdjacent(map->getWgtAboveLevel(d), 0, &counts[d][0]); 
      }
    }
    //calculate the overlap between this current event with the prob map.
    //got returned density and base_contribution
    O = map->overlap(sum,gps,AnitaPol::AnitaPol_t(pol),peak,true,S, &base_contribution, UCorrelator::ProbabilityMap::OVERLAP_SUM_SQRTS, !removed,0,&max_dens,&inv_two_pi_sqrt_det) / sqrt(S); 

    for (int level = 0; level < 10; level++){
      nclustered[level] = -1; 
      //for each level get the prob density threshold
      double thresh = UCorrelator::ProbabilityMap::dist2dens( map->getLevel(level), inv_two_pi_sqrt_det); 

      for (unsigned si = 0; si < max_dens.size(); si++){
        if (max_dens[si].second < thresh) continue; 

        if (counts[level][max_dens[si].first]){
           nclustered[level] = (int) counts[level][max_dens[si].first]; 
           break; 
        }
      }
    }

    if (mc) delete map; 
   
//    if (O < 0) O = -1; 
    printf("i=%d, run=%d, eventnumber=%d, F=%g, O=%g\n",i, run,eventNumber, F, O); 

    max_base_index = -1;
    base_sum = 0;
    max_base_p = 0;
    for (unsigned i = 0; i < base_contribution.size();i++){
      //for each base in base_contribution of this event
      // find the maximum prob density sum contribtion base. return that to max_base_index which is the base number.
      if (base_contribution[i].second > max_base_p){
        max_base_index = base_contribution[i].first;
        max_base_p = base_contribution[i].second; 
      }
      // base sum if suming of all ps from all bases for this events.
      base_sum += base_contribution[i].second; 
    }

    outputFile.cd(); 
    //all the branch variable are defined, so fill this event in output tree.
    outputTree->Fill(); 
  }

  outputFile.cd(); 
  outputTree->Write(); 
  return 0; 

}
