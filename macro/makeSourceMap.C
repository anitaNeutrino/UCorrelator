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
  p->collision_detection = false; 
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

int _makeSourceMap(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run = 50, int end_run =367, const char * outputFilePrefix = "source_maps/test")
{

  // Start getting the run / event numbers of events that pass our cut

  TChain c(treeName); 
  addRuns(c,start_run,end_run, thermalTreeFormat); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  int total_event_n = c.Draw("run:eventNumber:pol*5+peak:F",  TCut(TString::Format("(%s) && (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 
  printf("%d events pass selection\n", total_event_n); 


  UCorrelator::ProbabilityMap::Params * p = map_params(); 

  TFile f(TString::Format("%s%d_%d.root",outputFilePrefix,start_run, end_run), "RECREATE"); 
  // UCorrelator::ProbabilityMap *map_weighted = new UCorrelator::ProbabilityMap(p); 
  UCorrelator::ProbabilityMap map(p); 

  TTree tr("events","events"); 
  int run, ev, pol, peak,  nsegs ; 
  double S, F, p_ground, theta; 

  tr.Branch("event",&ev); 
  tr.Branch("run",&run); 
  tr.Branch("pol",&pol); 
  tr.Branch("peak",&peak); 
  tr.Branch("S",&S); 
  tr.Branch("F",&F);
  tr.Branch("p_ground",&p_ground);
  tr.Branch("theta",&theta);
  tr.Branch("nsegs",&nsegs);
  

  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  double last_integ = 0; 
  double last_integ_norm = 0; 
  // for (int i = 0; i < total_event_n; i+=1) 
  for (int i = 0; i < total_event_n; i+=1000) 
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
    // nsegs = map_weighted->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, S); 
    nsegs = map.add(p_ground, sum, gps, AnitaPol::AnitaPol_t(pol), peak, S);
    theta = -1*sum->peak[pol][peak].theta;
    // if(p_ground< 0.001){    
      printf("index = %d \t run = %d \t eventNumber = %d \t F = %g \t S = %g\t nsegs=%d \t p_ground = %g  theta= %g \n",i,run,ev,F,S,nsegs,p_ground, theta);
      // std::cout<< "\tsnr = "<< sum->deconvolved_filtered[pol][peak].snr << " longitude="<<sum->peak[pol][peak].longitude<<" latitude"<<sum->peak[pol][peak].latitude<< std::endl; 
    // }
    tr.Fill(); 

    // if (F > cutoff) 
    // {
    //   map.add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, 1); 
    // }
  }
  f.cd(); 
  map.Write("map_unweighted"); 
  // map_weighted->Write("map_weighted"); 
  tr.Write(); 
  return 0; 
}


int _evaluateSourceMap(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run = 50, int end_run =367, 
                       const char * sourceMapTree = "map_unweighted", const char * removed_events_file = "removed_events.txt", const char * outputFilePrefix = "source_maps_eval/",const char * sourcemapFilePrefix = "source_maps/test_") 
{
  
  TChain c(treeName); 
  addRuns(c,start_run,end_run, thermalTreeFormat); 
  TFile outputFile(TString::Format("%s%d_%d_%s.root",outputFilePrefix,start_run,end_run, treeName ),"RECREATE"); 
  //output file and output tree named overlap
  TTree outputTree("overlap","Overlap"); 
  double O=999,S,theta,base_sum,polangle,F,max_base_p; 
  int run,eventNumber,pol,peak,max_base_index;
  double wgt = 1; 
  double mcE = 0; 
  int nclustered[10]; 
  int removed = 0; 

  //branches in outputTree, the main purpose are to fill those branches.
  outputTree.Branch("O",&O); 
  outputTree.Branch("S",&S); 
  outputTree.Branch("run",&run); 
  outputTree.Branch("F",&F); 
  outputTree.Branch("eventNumber",&eventNumber); 
  outputTree.Branch("pol",&pol); 
  outputTree.Branch("theta",&theta); 
  outputTree.Branch("nclustered",&nclustered,"nclustered[10]/I"); 
  outputTree.Branch("polangle",&polangle); 
  outputTree.Branch("base_sum",&base_sum); //ps sum from all bases.
  outputTree.Branch("max_base_index",&max_base_index); // max ps base's id
  outputTree.Branch("max_base_p",&max_base_p); //max ps base's ps
  outputTree.Branch("weight",&wgt); // simulation event's weight 
  outputTree.Branch("mcE",&mcE); // simlation evetns's energy
  outputTree.Branch("peak",&peak); 
  outputTree.Branch("removed",&removed); // is this event removed or not.

  TFile sourceMapFile(TString::Format("%s%d_%d.root",sourcemapFilePrefix,start_run, end_run)); 
  UCorrelator::ProbabilityMap * map = (UCorrelator::ProbabilityMap*) sourceMapFile.Get(sourceMapTree); 
  UCorrelator::ProbabilityMap * source_map = map; //for mc 
  UCorrelator::ProbabilityMap::Params * map_pars = map_params(); 
  int total_event_n = c.Draw("run:eventNumber:pol*5+peak:F",  TCut(TString::Format("(%s) && (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 
  std::cout<< "total number of events: "<< total_event_n<< std::endl;
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
  for (int i = 0; i < total_event_n; i+=1){
    //loop through thermal chain
    run = c.GetV1()[i]; 
    if (run!= loaded_run)
    {
        if (sumfile) delete sumfile; 
        sumfile = new TFile(TString::Format(summaryFileFormat,run));        
        gROOT->cd(); 
        sumtree = (TTree*) sumfile->Get(treeName); 
        sumtree->SetBranchAddress("summary",&sum); 
        sumtree->SetBranchAddress("pat",&gps); 
        sumtree->BuildIndex("eventNumber"); 
        loaded_run =run; 
    }
    //each eventNumber in c, get S, pol, peak, F, eventNumber
    eventNumber = int(c.GetV2()[i]); 
    S = c.GetW()[i]; 
    pol = int(c.GetV3()[i]) / 5; 
    peak = int(c.GetV3()[i]) % 5;
    F = c.GetV4()[i]; 
    //get more variables from summary file
    sumtree->GetEntryWithIndex(eventNumber);  
    polangle = 90/TMath::Pi() * TMath::ATan2(sum->deconvolved_filtered[pol][peak].max_dU, sum->deconvolved_filtered[pol][peak].max_dQ); 
    theta = sum->peak[pol][peak].theta; 

    //this eventNumber is in the removed list.
    removed = (removed_events && removed_events->count(eventNumber));
    //fill the values for mc  
    wgt = sum->mc.weight; 
    mcE = sum->mc.energy; 


    if (F < cutoff) continue; 
    std::vector<std::pair<int,double> > base_contribution;
    //The maximun density of each bin.
    std::vector<std::pair<int,double> > max_dens;
    double inv_two_pi_sqrt_det; 

    // if (mc){
    //   // if mc, create new map, add this current event, combine with the initial source map.
    //   map = new UCorrelator::ProbabilityMap(map_pars); 
    //   map->add(sum,gps,AnitaPol::AnitaPol_t(pol),peak,S); 
    //   map->combineWith(*source_map); 
    //   for (int d = 0; d< 10; d++){
    //     map->groupAdjacent(map->getWgtAboveLevel(d), 0, &counts[d][0]); 
    //   }
    // }
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

    // if (mc) delete map; 
   
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
    outputTree.Fill(); 
  }

  outputFile.cd(); 
  outputTree.Write(); 
  return 0; 

}


void makeSourceMap(const char * treeName, bool evaluate = 1){
  int start_run,end_run;
  const char * sourceMapTree = "map_unweighted";
  if (!strcmp(treeName,"wais")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/wais/%d_max_30001_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/wais_%d-%d_max_30001_sinsub_10_3_ad_2.root";
    start_run = 120;
    end_run = 155;
    _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run);
    // _evaluateSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, sourceMapTree);

  }else if(!strcmp(treeName,"anita4")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/a4all/%d_max_30001_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/a4all_%d-%d_max_30001_sinsub_10_3_ad_2.root";
    start_run = 50;
    end_run = 367;
    // _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run);
    _evaluateSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, sourceMapTree);
  }else if(!strcmp(treeName,"simulation")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/simulated/%d_max_501_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/simulated_%d-%d_max_501_sinsub_10_3_ad_2.root";
    start_run = 1;
    end_run = 10;
    // _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run);
    _evaluateSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, sourceMapTree);
  }else{
    std::cout<< "wrong input treeName"<<std::endl;
  }
}




