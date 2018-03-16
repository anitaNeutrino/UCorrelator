#include "FFTtools.h" 
#include "AnitaConventions.h" 
#include "AnitaDataset.h" 

TCut cutString("theta<-3.5 && deconvImpulsivity>0.71");
// TCut cutString("theta<-3.5 && deconvImpulsivity>0.71");
// const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && theta < 40 )";
// const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && isMostImpulsive && !payloadBlast && MaxPeak < 1000 && theta < 40 && ( (HPolTrigger && iteration < 5) || (VPolTrigger && iteration > 4))  && !isPulser  )";

UCorrelator::ProbabilityMap::Params * map_params()
{
  // pixel x, pixel y , max meter x, max meter y
  StereographicGrid * g= new StereographicGrid(1000,1000,2000000,2000000); 

  TF1 * f_dtheta = new TF1("ftheta", "[0]/x^[1] + [2]", 1, 100);
  TF1 * f_dphi = new TF1("fphi", "[0]/x^[1] + [2]", 1, 100);
  //anita4 fit from wais.
  f_dtheta->SetParameter(0, 4.714); 
  f_dtheta->SetParameter(1, 1.211); 
  f_dtheta->SetParameter(2, 0.06639); 
  f_dphi->SetParameter(0, 40.96); 
  f_dphi->SetParameter(1, 1.61);
  f_dphi->SetParameter(2, 0.1766);
  //anita3
  // f_dtheta->SetParameter(0, 0.3936); 
  // f_dtheta->SetParameter(1, 0.2102); 
  // f_dphi->SetParameter(0, 1.065); 
  // f_dphi->SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel * snrResolutionModel = new UCorrelator::PointingResolutionParSNRModel (*f_dtheta, *f_dphi, true, true); // the last false is not to use cos_theta_scale
  // UCorrelator::PointingResolutionModelPlusHeadingError * resolutionModel = new UCorrelator::PointingResolutionModelPlusHeadingError(20, snrResolutionModel); 

  Refraction::SphRay * ref = new Refraction::SphRay; 

  UCorrelator::ProbabilityMap::Params *p = new UCorrelator::ProbabilityMap::Params; 
  p->refract = ref; 
  // p->refract = 0; 
  p->seg = g; 
  p->point = snrResolutionModel; 
  p->collision_detection = false; 
  p->verbosity = 0; // verbosity level for output info.
  p->maximum_distance = 2.51;
  // p->min_p_on_continent = 0;
 

  return p; 

}
void clearTheObject(UCorrelator::ProbabilityMap::Params * p){
  delete p->seg; 
  delete p->point;
  delete p;
}


int combineSourceMaps(const char * dir, const char * output)
{

  TSystemDirectory d(dir,dir); 
  TList * files = d.GetListOfFiles(); 

  TSystemFile * file;
  TString fname; 
  TIter next(files); 

  UCorrelator::ProbabilityMap * map_u = 0;
  TFile outf(output,"RECREATE");
  //loop through the files in the dir, add map together and put a new map in output folder. 
  while ((file=(TSystemFile*) next()))
  {
    fname = file->GetName(); 
    fname = TString(dir) + TString("/") + fname; 
    if (fname.EndsWith(".root"))
    {
      TFile * f = new TFile(fname); 
      printf("Considering %s\n", fname.Data());       
      UCorrelator::ProbabilityMap * m_u = (UCorrelator::ProbabilityMap*) f->Get("map_unweighted");
      if (!map_u)
      {
        map_u = m_u; 
      }
      else
      {
        map_u->combineWith(*m_u); 
        delete m_u; 
      }
    }
  }

  outf.cd(); 
  map_u->Write("map_unweighted"); 

  return 0;
}
int _trendOfSinglets(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run, int end_run, const char * filePrefix, int mod){
  
  UCorrelator::ProbabilityMap * map = 0;
  UCorrelator::ProbabilityMap * removeMap = 0;
  UCorrelator::ProbabilityMap * combineMap = 0;
  std::vector<UCorrelator::ProbabilityMap *> maps;
  for(int i = 0; i < mod; i++){
    std::cout << "loaded all the map "<< i << std::endl;
    TFile * tempFile = TFile::Open(TString::Format("source_maps_withSNRcut/%s/%smod%d_remainder%d_%d_%d.root",treeName,filePrefix,mod, i ,start_run, end_run)); 
    maps.push_back((UCorrelator::ProbabilityMap*) tempFile->Get("map_unweighted"));
    tempFile->Close();
    // delete tempFile;
  }
  std::cout << "loaded all the maps once for all!!!!!!!!!" << std::endl;
  // return 0;
  TFile* removeFile;
  TFile* combineFile;

  TFile outputFile(TString::Format("source_maps_withSNRcut/%s/trendOfSinglets%smod%d_%d_%d.root",treeName,filePrefix,mod,start_run,end_run),"RECREATE"); 
  //output file and output tree named Overlap
  TTree outputTree("trend","trend"); 
  double percentOfData;
  double N_singlets;  
  double N_singlets_nearbase;  
  double N_singlets_notnearbase;  
  //branches in outputTree, the main purpose are to fill those branches.
  outputTree.Branch("percentOfData",&percentOfData); 
  outputTree.Branch("N_singlets",&N_singlets); 
  outputTree.Branch("N_singlets_nearbase",&N_singlets_nearbase); 
  outputTree.Branch("N_singlets_notnearbase",&N_singlets_notnearbase); 
      
  for(int l =0; l<mod; l++){
    int totalSinglets = 0;
    int totalSingletsNearBase = 0;
    int totalSingletsNotNearBase = 0;
    for (int i = 0; i < mod; i++){
      if(i == 0 and l == 0){
        //start 
        std::cout <<"add file index "<< 0 << std::endl;
        combineFile = new TFile(TString::Format("source_maps_withSNRcut/%s/%smod%d_remainder%d_%d_%d.root",treeName,filePrefix,mod, 0 ,start_run, end_run)); 
        map = (UCorrelator::ProbabilityMap*) combineFile->Get("map_unweighted");
      }else{ 
        if (i == 0){
          std::cout <<"add file index "<< (i + l - 1)%mod << std::endl;
          combineMap = maps[(i + l - 1)%mod];
          map->combineWith(*combineMap);
        }
        //otherwise, just need to remove the i-1 file and add (i+l)%mod file.
        std::cout <<"remove file index "<< (i - 1 + mod)%mod<< std::endl;
        removeMap = maps[(i - 1 + mod)%mod];
        map->removeWith(*removeMap);

        std::cout <<"add file index "<< (i + l)%mod << std::endl;
        combineMap = maps[(i + l)%mod];
        map->combineWith(*combineMap);
      } 

      std::pair<int, int> results = map->showClusters(0,0);
      totalSingletsNearBase += results.first;
      totalSingletsNotNearBase += results.second;
      totalSinglets += results.first + results.second;
      std::cout << "N_singlets "<< results.first + results.second << " totalSinglets "<< totalSinglets << std::endl;
      std::cout<<"l "<<l <<" i " << i<< " events "<< map->getProbSumsIntegral(true) <<  std::endl;
    }
    std::cout <<"====== avg nSinglets" << totalSinglets/float(mod) <<" totalSinglets "<<totalSinglets<< " mod "<< mod << " numberOfFiles="<< l + 1 << std::endl;
    N_singlets_nearbase = totalSingletsNearBase/float(mod);  
    N_singlets_notnearbase = totalSingletsNotNearBase/float(mod);  
    N_singlets = totalSinglets/float(mod);  
    percentOfData = (l + 1)/float(mod);
    outputTree.Fill();
  }
  outputFile.Write();
  delete map;
  return 0;
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

int _makeSourceMap(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run = 50, int end_run =367, const char * filePrefix = "_2.51sigma_1pc_", int mod=1, int mod_remainder=0)
{

  // Start getting the run / event numbers of events that pass our cuts

  TChain c(treeName); 
  addRuns(c,start_run,end_run, thermalTreeFormat); 
  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  int total_event_n = c.Draw("run:eventNumber:pol*5+peak:deconvImpulsivity",  cutString, "goff" ); 
  printf("%d events pass selection\n", total_event_n); 


  UCorrelator::ProbabilityMap::Params * p = map_params(); 

  TFile f(TString::Format("source_maps_withSNRcut/%s/%smod%d_remainder%d_%d_%d.root",treeName,filePrefix,mod, mod_remainder,start_run, end_run), "RECREATE"); 
  // UCorrelator::ProbabilityMap *map_weighted = new UCorrelator::ProbabilityMap(p); 
  UCorrelator::ProbabilityMap map(p); 

  TTree tr("events","events"); 
  int run, ev, pol, peak,  nsegs ; 
  double S, deconvImpulsivity, p_ground, theta,snr; 
  tr.Branch("event",&ev); 
  tr.Branch("run",&run); 
  tr.Branch("pol",&pol); 
  tr.Branch("peak",&peak); 
  tr.Branch("S",&S); 
  tr.Branch("deconvImpulsivity",&deconvImpulsivity);
  tr.Branch("p_ground",&p_ground);
  tr.Branch("theta",&theta);
  tr.Branch("nsegs",&nsegs);
  tr.Branch("snr",&snr);
  

  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  double last_integ = 0; 
  double last_integ_norm = 0; 
  // for (int i = mod_remainder; i < total_event_n; i+=mod) 
  for(int i =0; i< total_event_n; i++){
    if (TString::Hash(&i, sizeof(i))%mod != mod_remainder){
      continue;
    }
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
    // if(pol == 0){
    //   continue;
    // } 
    peak = int(c.GetV3()[i]) % 5; 
    deconvImpulsivity = c.GetV4()[i]; 
    // nsegs = map_weighted->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, S); 
    nsegs = map.add(p_ground, sum, gps, AnitaPol::AnitaPol_t(pol), peak, S);
    theta = -1*sum->peak[pol][peak].theta;
    snr = sum->peak[pol][peak].snr;
    // if(p_ground< 0.1){    
      printf("index = %d \t run = %d \t eventNumber = %d \t deconvImpulsivity = %g \t S = %g\t nsegs=%d \t p_ground = %g  theta= %g \n",i,run,ev,deconvImpulsivity,S,nsegs,p_ground, theta);
      // std::cout<< "\tsnr = "<< sum->deconvolved_filtered[pol][peak].snr << " longitude="<<sum->peak[pol][peak].longitude<<" latitude"<<sum->peak[pol][peak].latitude<< std::endl; 
    // }
    tr.Fill(); 
  }
  f.cd(); 
  map.Write("map_unweighted"); 
  // map_weighted->Write("map_weighted"); 
  tr.Write();
  clearTheObject(p);
  return 0; 
}


int _evaluateSourceMap(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run, int end_run, const char * filePrefix, int mod, int mod_remainder) 
{
  const char * sourceMapTree = "map_unweighted";
  TChain c(treeName); 
  addRuns(c,start_run,end_run, thermalTreeFormat); 
  TFile outputFile(TString::Format("source_maps_eval/%s/%smod%d_remainder%d_%d_%d.root",treeName,filePrefix,mod,mod_remainder,start_run,end_run),"RECREATE"); 
  //output file and output tree named overlap
  TTree outputTree("overlap","Overlap"); 
  double O=999,S,theta,base_sum,polangle,deconvImpulsivity,max_base_p,snr,longitude,latitude; 
  int run,eventNumber,pol,peak,max_base_index;
  double wgt = 1; 
  double mcE = 0; 
  int nclustered[10]; 
  int removed = 0; 

  //branches in outputTree, the main purpose are to fill those branches.
  outputTree.Branch("O",&O); 
  outputTree.Branch("S",&S); 
  outputTree.Branch("run",&run); 
  outputTree.Branch("deconvImpulsivity",&deconvImpulsivity); 
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
  outputTree.Branch("snr",&snr); 
  outputTree.Branch("longitude",&longitude); 
  outputTree.Branch("latitude",&latitude); 
  
  TFile sourceMapFile(TString::Format("source_maps_withSNRcut/%s/%smod%d_remainder%d_%d_%d.root",treeName,filePrefix,mod, mod_remainder,start_run, end_run)); 
  UCorrelator::ProbabilityMap * map = (UCorrelator::ProbabilityMap*) sourceMapFile.Get(sourceMapTree); 
  UCorrelator::ProbabilityMap * source_map = map; //for mc 
  UCorrelator::ProbabilityMap::Params * map_pars = map_params(); 
  int total_event_n = c.Draw("run:eventNumber:pol*5+peak:deconvImpulsivity",  cutString, "goff" ); 
  std::cout<< "total number of events: "<< total_event_n<< std::endl;
  std::vector<std::vector<double> > counts(10, std::vector<double> (map->segmentationScheme()->NSegments())); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 
  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  // for (int i = mod_remainder; i < total_event_n; i+=mod){
  for(int i =0; i< total_event_n; i++){
    //loop through thermal chain
    if (TString::Hash(&i, sizeof(i))%mod != mod_remainder){
      continue;
    }
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
    //each eventNumber in c, get S, pol, peak, deconvImpulsivity, eventNumber
    eventNumber = int(c.GetV2()[i]); 
    S = c.GetW()[i]; 
    pol = int(c.GetV3()[i]) / 5; 
    // if(pol == 0){
    //   continue;
    // } 
    peak = int(c.GetV3()[i]) % 5;
    deconvImpulsivity = c.GetV4()[i]; 
    //get more variables from summary file
    sumtree->GetEntryWithIndex(eventNumber);  
    polangle = 90/TMath::Pi() * TMath::ATan2(sum->deconvolved_filtered[pol][peak].max_dU, sum->deconvolved_filtered[pol][peak].max_dQ); 
    theta = -1*sum->peak[pol][peak].theta; 
    snr = sum->peak[pol][peak].snr; 
    longitude = sum->peak[pol][peak].longitude;
    latitude = sum->peak[pol][peak].latitude;
    //fill the values for mc  
    wgt = sum->mc.weight; 
    mcE = sum->mc.energy; 
    std::vector<std::pair<int,double> > base_contribution;
    //The maximun density of each bin.
    std::vector<std::pair<int,double> > max_dens;
    double inv_two_pi_sqrt_det; 

    //calculate the overlap between this current event with the prob map.
    //got returned density and base_contribution
    O = map->overlap(sum,gps,AnitaPol::AnitaPol_t(pol),peak,true,S, &base_contribution, UCorrelator::ProbabilityMap::OVERLAP_SUM_SQRTS, !removed,0,&max_dens,&inv_two_pi_sqrt_det) / sqrt(S); 
    // if (mc) delete map; 
   
//    if (O < 0) O = -1; 
      printf("i=%d, run=%d, eventnumber=%d, deconvImpulsivity=%g, O=%g\n",i, run,eventNumber, deconvImpulsivity, O); 

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
  clearTheObject(map_pars);
  return 0; 

}


void makeSourceMap(const char * treeName, bool evaluate = 1){
  int start_run,end_run;
  const char * sourceMapTree = "map_unweighted";
  if (!strcmp(treeName,"wais")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/wais/%d_max_30001_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/wais_%d-%d_max_30001_sinsub_10_3_ad_2.root";
    const char * filePrefix = "_2.51sigma_30001_";
    int mod = 1;
    int mod_remainder = 0;
    start_run = 120;
    end_run = 155;
    _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
    // _evaluateSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);

  }else if(!strcmp(treeName,"anita4")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/a4all/%d_max_30002_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/a4all_%d-%d_max_30002_sinsub_10_3_ad_2.root";
    const char * filePrefix = "_2.51sigma_30002_";
    int mod = 1;
    int mod_remainder = 0;
    start_run = 40;
    end_run = 367;
    for (mod_remainder= 0; mod_remainder<mod; mod_remainder++){
    // for (mod_remainder= 1; mod_remainder<2; mod_remainder++){
       _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
       // _evaluateSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
    }
    // _trendOfSinglets(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod);

  }else if(!strcmp(treeName,"simulation")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/simulated/%d_max_1001_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/simulated_%d-%d_max_1001_sinsub_10_3_ad_2.root";
    const char * filePrefix = "_2.51sigma_1001_";
    int mod = 1;
    int mod_remainder = 0;
    start_run = 1;
    end_run = 10;
    _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
    // _evaluateSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
  }else{
    std::cout<< "wrong input treeName"<<std::endl;
  }
}




