#include "FFTtools.h" 
#include "AnitaConventions.h" 
#include "AnitaDataset.h" 
 
// TCut cutString("");
TCut cutString("theta<-3.5 && impulsivity>0.71");
double xSigma = 5;
const char * xString = "5";
// TCut cutString("theta<-5.8 && impulsivity<0.60"); // select for thermal events
// TCut cutString("theta<-3.5 && impulsivity>0.71");
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
  UCorrelator::PointingResolutionParSNRModel * snrResolutionModel = new UCorrelator::PointingResolutionParSNRModel (*f_dtheta, *f_dphi, true, false); // the last false is not to use cos_theta_scale
  // UCorrelator::PointingResolutionModelPlusHeadingError * resolutionModel = new UCorrelator::PointingResolutionModelPlusHeadingError(20, snrResolutionModel); 

  Refraction::SphRay * ref = new Refraction::SphRay; 

  UCorrelator::ProbabilityMap::Params *p = new UCorrelator::ProbabilityMap::Params; 
  p->refract = ref; 
  // p->refract = 0; 
  p->seg = g; 
  p->point = snrResolutionModel; 
  p->collision_detection = false; 
  p->verbosity = 0; // verbosity level for output info.
  p->maximum_distance = xSigma;
  p->max_dphi = 5;
  p->max_dtheta = 5;
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

  UCorrelator::ProbabilityMap * map = 0;
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
      UCorrelator::ProbabilityMap * m_u = (UCorrelator::ProbabilityMap*) f->Get("maps");
      if (!map)
      {
        map = m_u; 
      }
      else
      {
        map->combineWith(*m_u); 
        delete m_u; 
      }
    }
  }

  outf.cd(); 
  map->Write("maps"); 

  return 0;
}
int _trendOfSinglets(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run, int end_run, const char * filePrefix, int mod){
  
  UCorrelator::ProbabilityMap * map = 0;
  UCorrelator::ProbabilityMap * removeMap = 0;
  UCorrelator::ProbabilityMap * combineMap = 0;
  std::vector<UCorrelator::ProbabilityMap *> maps;
  for(int i = 0; i < mod; i++){
    std::cout << "loaded the map "<< i << std::endl;
    TFile * tempFile = TFile::Open(TString::Format("source_maps/%s/_%s%smod%d_remainder%d_%d_%d.root",treeName,xString,filePrefix,mod, i ,start_run, end_run)); 
    maps.push_back((UCorrelator::ProbabilityMap*) tempFile->Get("maps"));
    tempFile->Close();
    // delete tempFile;
  }
  std::cout << "loaded all the maps once for all!!!!!!!!!" << std::endl;
  // return 0;
  TFile* removeFile;
  TFile* combineFile;

  TFile outputFile(TString::Format("source_maps/%s/trendOfSinglets_%s%smod%d_%d_%d.root",treeName,xString,filePrefix,mod,start_run,end_run),"RECREATE"); 
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
        combineFile = new TFile(TString::Format("source_maps/%s/_%s%smod%d_remainder%d_%d_%d.root",treeName,xString,filePrefix,mod, 0 ,start_run, end_run)); 
        map = (UCorrelator::ProbabilityMap*) combineFile->Get("maps");
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

int _makeSourceMap(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run = 50, int end_run =367, const char * filePrefix = "", int mod=1, int mod_remainder=0, int weighted = 0)
{

  // Start getting the run / event numbers of events that pass our cuts

  TChain c(treeName); 
  addRuns(c,start_run,end_run, thermalTreeFormat); 
  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  int total_event_n = c.Draw("run:eventNumber:pol*5+peak:impulsivity:linearPolFrac:linearPolAngle:powerH:powerV",  cutString, "goff" ); 
  printf("%d events pass selection\n", total_event_n); 


  UCorrelator::ProbabilityMap::Params * p = map_params(); 
  TFile f(TString::Format("source_maps/%s/_%s%smod%d_remainder%d_%d_%d.root",weighted ? "weightedMC":treeName,xString,filePrefix,mod, mod_remainder,start_run, end_run), "RECREATE"); 
 
  // UCorrelator::ProbabilityMap *map_weighted = new UCorrelator::ProbabilityMap(p); 
  UCorrelator::ProbabilityMap map(p); 

  TTree tr("events","events"); 
  int run, ev, pol, peak,  nsegs,  NOverlapedBases ; 
  double S, impulsivity,impulsivityH,impulsivityV,powerH,powerV,linearPolFrac,linearPolAngle, p_ground, theta,phi,snr,longitude,latitude; 
  tr.Branch("event",&ev); 
  tr.Branch("run",&run); 
  tr.Branch("pol",&pol); 
  tr.Branch("peak",&peak); 
  tr.Branch("S",&S); 
  tr.Branch("impulsivity",&impulsivity);
  tr.Branch("impulsivityH",&impulsivityH);
  tr.Branch("impulsivityV",&impulsivityV);
  tr.Branch("powerH",&powerH);
  tr.Branch("powerV",&powerV);
  tr.Branch("linearPolFrac",&linearPolFrac);
  tr.Branch("linearPolAngle",&linearPolAngle);
  tr.Branch("NOverlapedBases",&NOverlapedBases);
  tr.Branch("p_ground",&p_ground);
  tr.Branch("theta",&theta);
  tr.Branch("phi",&phi);
  tr.Branch("nsegs",&nsegs);
  tr.Branch("snr",&snr);
  tr.Branch("longitude",&longitude); 
  tr.Branch("latitude",&latitude); 
  tr.Branch("gps",&gps); 
  

  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  double last_integ = 0; 
  double last_integ_norm = 0; 
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
    sumtree->GetEntryWithIndex(ev); 
    pol = int(c.GetV3()[i]) / 5;
    // if(pol == 0){
    //   continue;
    // } 
    peak = int(c.GetV3()[i]) % 5; 
    impulsivity = c.GetV4()[i]; 
    // nsegs = map_weighted->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, S); 
    nsegs = map.add(NOverlapedBases, p_ground, sum, gps, AnitaPol::AnitaPol_t(pol), peak, S);
    if(pol == 0){
      impulsivityH = sum->deconvolved_filtered[0][peak].impulsivityMeasure;
      impulsivityV = sum->deconvolved_filtered[1][0].impulsivityMeasure;
    }else{
      impulsivityH = sum->deconvolved_filtered[0][0].impulsivityMeasure;
      impulsivityV = sum->deconvolved_filtered[1][peak].impulsivityMeasure;
    }
    linearPolFrac = c.GetVal(4)[i];
    linearPolAngle = c.GetVal(5)[i];
    powerH = c.GetVal(6)[i];
    powerV = c.GetVal(7)[i];
    theta = -1*sum->peak[pol][peak].theta;
    phi = sum->peak[pol][peak].phi;
    snr = sum->peak[pol][peak].snr;
    longitude = sum->peak[pol][peak].longitude;
    latitude = sum->peak[pol][peak].latitude;
    if(weighted){
      S = sum->mc.weight;
    }else{
      S = 1;
    }
    // if(p_ground< 0.1){    
      printf("index= %d \trun = %d \tevn= %d \timpulsivity= %g \tS= %g\t nsegs= %d \tp_ground= %g \ttheta= %g \tsnr=%g \n",i,run,ev,impulsivity,S,nsegs,p_ground, theta, snr);
      // std::cout<< "\tsnr = "<< sum->deconvolved_filtered[pol][peak].snr << " longitude="<<sum->peak[pol][peak].longitude<<" latitude"<<sum->peak[pol][peak].latitude<< std::endl; 
    // }
    tr.Fill(); 
  }
  //doClustering
  map.doClustering(); 
  //each event find it belongs to which cluster.

  double indexOfCluster,sizeOfCluster;
  TBranch *branchIndexOfCluster = tr.Branch("indexOfCluster",&indexOfCluster);
  TBranch *branchSizeOfCluster = tr.Branch("sizeOfCluster",&sizeOfCluster);
  tr.SetBranchAddress("theta",&theta);
  tr.SetBranchAddress("phi",&phi);
  tr.SetBranchAddress("gps",&gps);
  tr.SetBranchAddress("nsegs",&nsegs);
  for(int i =0; i< tr.GetEntries(); i++){
    tr.GetEntry(i);
    if(nsegs!=0){
       map.evaluateEvent(indexOfCluster, sizeOfCluster, -1*theta,phi, gps);
    }else{
      indexOfCluster = -1;
      sizeOfCluster = -1;
    }
    printf("%g ",indexOfCluster);
    // indexOfCluster = -1 :  the event does not project to ground with certern sigma, or nsegs == 0  
    branchIndexOfCluster->Fill();
    branchSizeOfCluster->Fill();
  }
  f.cd(); 
  map.Write("maps"); 
  tr.Write();
  clearTheObject(p);
  return 0; 
}


int _evaluateSourceMap(const char * treeName, const char* summaryFileFormat, const char* thermalTreeFormat, int start_run, int end_run, const char * filePrefix, int mod, int mod_remainder) 
{
std::cout<<"hello 0"<<std::endl;

  TFile sourceMapFile(TString::Format("source_maps/%s/_%s%smod%d_remainder%d_%d_%d.root",treeName,xString,filePrefix,mod, mod_remainder,start_run, end_run)); 
  TTree * events =(TTree *) sourceMapFile.Get("events"); 
  int total_event_n = events->GetEntries();
  std::cout<< "total number of events: "<< total_event_n<< std::endl;
  int run, event, pol, nsegs, NOverlapedBases;
  double impulsivity, impulsivityV, impulsivityH, powerV, powerH, linearPolFrac, linearPolAngle, indexOfCluster, longitude, latitude;
  events->SetBranchAddress("run",&run);
  events->SetBranchAddress("event",&event);
  events->SetBranchAddress("pol",&pol);
  events->SetBranchAddress("nsegs",&nsegs);
  events->SetBranchAddress("indexOfCluster",&indexOfCluster);
  events->SetBranchAddress("NOverlapedBases",&NOverlapedBases);
  events->SetBranchAddress("impulsivity",&impulsivity);
  events->SetBranchAddress("impulsivityV",&impulsivityV);
  events->SetBranchAddress("impulsivityH",&impulsivityH);
  events->SetBranchAddress("powerV",&powerV);
  events->SetBranchAddress("powerH",&powerH);
  events->SetBranchAddress("linearPolFrac",&linearPolFrac);
  events->SetBranchAddress("linearPolAngle",&linearPolAngle);
  events->SetBranchAddress("longitude",&longitude);
  events->SetBranchAddress("latitude",&latitude);
  int count_n[39000]={0}, count_base[39000]={0}, count_noBase[39000]={0}, count_H[39000]={0}, count_Mix[39000]={0}, count_V[39000]={0};
  int count_baseH[39000]={0}, count_baseMix[39000]={0}, count_baseV[39000]={0}, count_noBaseH[39000]={0}, count_noBaseMix[39000]={0}, count_noBaseV[39000]={0};
  int first_pol[39000] = {0};
  double sum_linearPolFrac[39000] = {0};
  double sum_linearPolAngle[39000] = {0};
  double sum_powerH[39000] = {0};
  double sum_powerV[39000] = {0};
  double first_longitude[39000] = {0};
  double first_latitude[39000] = {0};
  for(int i =0; i< total_event_n; i++){
    events->GetEntry(i);
    int j = round(indexOfCluster);
    if (nsegs == 0 or j <= 0){
      continue;
    }
    count_n[j]++ ;
    sum_linearPolFrac[j] += linearPolFrac;
    sum_linearPolAngle[j] += linearPolAngle;
    sum_powerH[j] += powerH;
    sum_powerV[j] += powerV;

    if(first_longitude[j] == 0){
      first_longitude[j] = longitude;
    }
    if(first_latitude[j] == 0){
      first_latitude[j] = latitude;
    }
    if(first_pol[j] == 0){
      first_pol[j] = pol;
    }
    // this event overlap with some base
    if(NOverlapedBases != 0){
      count_base[j] ++;
      if(pol==0){
        count_H[j]++;
        count_baseH[j]++;
      }else{
        count_V[j]++;
        count_baseV[j]++;
      }
      // if(impulsivityV - impulsivityH > 0.2){
      //   count_V[j]++;
      //   count_baseV[j]++;
      // }else if(impulsivityV - impulsivityH < -0.2){
      //   count_H[j]++;
      //   count_baseH[j]++;
      // }else{
      //   count_Mix[j]++;
      //   count_baseMix[j]++;
      // }
    }else{
      // this event overlap with no base
      count_noBase[j] ++;
      if(pol==0){
        count_H[j]++;
        count_noBaseH[j]++;
      }else{
        count_V[j]++;
        count_noBaseV[j]++;
      }
      // if(impulsivityV - impulsivityH > 0.2){
      //   count_V[j]++;
      //   count_noBaseV[j]++;
      // }else if(impulsivityV - impulsivityH < -0.2){
      //   count_H[j]++;
      //   count_noBaseH[j]++;
      // }else{
      //   count_Mix[j]++;
      //   count_noBaseMix[j]++;
      // }
    }
  }
  // prepare the output file
  TFile outputFile(TString::Format("cluster/%s/_%s%smod%d_remainder%d_%d_%d.root",treeName,xString,filePrefix,mod,mod_remainder,start_run,end_run),"RECREATE"); 
  //output file and output tree named cluster
  TTree outputTree("cluster","cluster"); 
  int index, n, base, noBase, H, Mix, V, baseH, baseMix, baseV, noBaseH, noBaseMix, noBaseV, _pol;
  double avgLinearPolFrac, avgLinearPolAngle,_powerH, _powerV, _longitude, _latitude;
  //branches in outputTree, the main purpose are to fill those branches.
  outputTree.Branch("index",&index); 
  outputTree.Branch("pol",&_pol); 
  outputTree.Branch("n",&n); 
  outputTree.Branch("base",&base); 
  outputTree.Branch("noBase",&noBase); 
  outputTree.Branch("H",&H); 
  outputTree.Branch("Mix",&Mix); 
  outputTree.Branch("V",&V); 
  outputTree.Branch("baseH",&baseH); 
  outputTree.Branch("baseMix",&baseMix); 
  outputTree.Branch("baseV",&baseV); 
  outputTree.Branch("noBaseH",&noBaseH); 
  outputTree.Branch("noBaseMix",&noBaseMix); 
  outputTree.Branch("noBaseV",&noBaseV); 
  outputTree.Branch("avgLinearPolFrac",&avgLinearPolFrac); 
  outputTree.Branch("avgLinearPolAngle",&avgLinearPolAngle);
  outputTree.Branch("powerH",&_powerH); 
  outputTree.Branch("powerV",&_powerV);  
  outputTree.Branch("longitude",&_longitude); 
  outputTree.Branch("latitude",&_latitude); 

  outputFile.cd(); 
  for(int j = 0; j< 39000; j++){
  // std::cout<<j<<" "<< count_n[j]<<std::endl;
    if(count_n[j]==0){
      continue;
    }
    index=j;
    n=count_n[j];
    _pol=first_pol[j];
    base=count_base[j];
    noBase=count_noBase[j];
    H=count_H[j];
    Mix=count_Mix[j];
    V=count_V[j];
    baseH=count_baseH[j];
    baseMix=count_baseMix[j];
    baseV=count_baseV[j];
    noBaseH=count_noBaseH[j];
    noBaseMix=count_noBaseMix[j];
    noBaseV=count_noBaseV[j];
    avgLinearPolFrac=sum_linearPolFrac[j]/double(n);
    avgLinearPolAngle=sum_linearPolAngle[j]/double(n);
    _powerH=sum_powerH[j]/double(n);
    _powerV=sum_powerV[j]/double(n);
    _longitude=first_longitude[j];
    _latitude=first_latitude[j];
    outputTree.Fill(); 
  }

    //all the branch variable are defined, so fill this event in output tree.
  outputFile.cd(); 
  outputTree.Write();
  return 0; 

}


void _makeMCmapAndEvaluateEfficiency(){
  char const * filePrefixs[15] = {"0.0002","1.0","1.8","2.0","2.2","2.5","3.5","4.5","5.5","6","6.5","7","7.5","8","10.5"};
  TFile * sourceMapFile;
  TFile mcMapFile(TString::Format("source_maps/weightedMC/weightedMCmod1_remainder0_1_500.root")); 
  UCorrelator::ProbabilityMap * mcMaps = (UCorrelator::ProbabilityMap*) mcMapFile.Get("maps");
  
  UCorrelator::ProbabilityMap * sourceMaps;
  for (int i=0; i<15; i++){
    sourceMapFile = new TFile(TString::Format("source_maps/anita4/_%ssigma_30002_mod1_remainder0_41_367.root",filePrefixs[i])); 
    sourceMaps = (UCorrelator::ProbabilityMap*) sourceMapFile->Get("maps");
    double nMasked=0, nNotMasked=0;
    mcMaps->maskingWithMap(nMasked, nNotMasked, *sourceMaps);
    std::cout<< "i="<<i<<"\t nMasked = "<<nMasked<< " \tnNotMasked = "<<nNotMasked<<"\t "<<filePrefixs[i]<< " \t"<< nNotMasked/(nMasked+nNotMasked)<<std::endl;
    delete sourceMapFile;
  }
  
}


void makeSourceMap(const char * treeName, bool evaluate = 1){
  int start_run,end_run;
  const char * sourceMapTree = "maps";
  if (!strcmp(treeName,"wais")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/wais/%d_max_30001_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/wais_%d-%d_max_30001_sinsub_10_3_ad_2.root";
    const char * filePrefix = "sigma_30001_";
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
    const char * filePrefix = "sigma_30002_";
    int mod = 1;
    int mod_remainder = 0;
    start_run = 41;
    end_run = 367;
    for (mod_remainder= 0; mod_remainder<mod; mod_remainder++){
    // for (mod_remainder= 1; mod_remainder<2; mod_remainder++){
       _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
       _evaluateSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
       // _makeMCmapAndEvaluateEfficiency(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
    }
    // _trendOfSinglets(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod);

  }else if(!strcmp(treeName,"simulation")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/simulated/%d_max_1000_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/simulated_%d-%d_max_1000_sinsub_10_3_ad_2.root";
    const char * filePrefix = "sigma_1000_";
    int mod = 170;
    int mod_remainder = 0;
    start_run = 1;
    end_run = 1;
    for (mod_remainder= 0; mod_remainder<mod; mod_remainder+=2){
    // for (mod_remainder= 0; mod_remainder<1; mod_remainder++){
      _makeSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
      _evaluateSourceMap(treeName, summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
    }
  }else if(!strcmp(treeName,"weightedMC")){
    std::cout<<"makeSourceMap: "<< treeName <<std::endl;
    const char* summaryFileFormat = "/Volumes/SDCard/data/simulated/%d_max_1000_sinsub_10_3_ad_2.root";
    const char* thermalTreeFormat = "thermalTrees/simulated_%d-%d_max_1000_sinsub_10_3_ad_2.root";
    const char * filePrefix = "sigma_1000_";
    int mod = 1;
    int mod_remainder = 0;
    start_run = 1;
    end_run = 500;
    // for (mod_remainder= 0; mod_remainder<mod; mod_remainder+=2){
    for (mod_remainder= 0; mod_remainder<1; mod_remainder++){
      // _makeSourceMap("simulation", summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder,1);
      // _evaluateSourceMap("simulation", summaryFileFormat, thermalTreeFormat, start_run, end_run, filePrefix, mod, mod_remainder);
    }
    _makeMCmapAndEvaluateEfficiency();

  }else{
    std::cout<< "wrong input treeName"<<std::endl;
  }
}




