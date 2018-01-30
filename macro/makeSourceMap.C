#include "FFTtools.h" 
#include "AnitaConventions.h" 
#include "AnitaDataset.h" 



double cutoff = 3; 
const char * weight = "F > 3";
// const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && theta < 40 )";
// const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && isMostImpulsive && !payloadBlast && MaxPeak < 1000 && theta < 40 && ( (HPolTrigger && iteration < 5) || (VPolTrigger && iteration > 4))  && !isPulser  )";



int combineSourceMaps(const char * dir, const char * output)
{

  TSystemDirectory d(dir,dir); 
  TList * files = d.GetListOfFiles(); 

  TSystemFile * file;
  TString fname; 
  TIter next(files); 

  UCorrelator::ProbabilityMap * map_w = 0;
  UCorrelator::ProbabilityMap * map_u = 0;

  TFile outf(output,"RECREATE"); 
  while ((file=(TSystemFile*) next()))
  {
    fname = file->GetName(); 
    fname = TString(dir) + TString("/") + fname; 
    if (fname.EndsWith(".root"))
    {
      TFile * f = new TFile(fname); 
      printf("Considering %s\n", fname.Data()); 
      
      UCorrelator::ProbabilityMap * m_w = (UCorrelator::ProbabilityMap*) f->Get("map_weighted");
      UCorrelator::ProbabilityMap * m_u = (UCorrelator::ProbabilityMap*) f->Get("map_unweighted");
      if (!m_w || !m_u) continue; 

      if (!map_w)
      {
        map_w = m_w; 
      }
      else
      {
        map_w->combineWith(*m_w); 
        delete m_w; 
      }

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
  map_w->Write("map_weighted"); 
  map_u->Write("map_unweighted"); 

  return 0;
}

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



void addRuns(TChain & c, int start_run, int end_run)
{
  for (int i = start_run; i <= end_run; i+=40) 
  {
    if (start_run < i + 40 && end_run >= i)
    {
      // TString adding = TString::Format("thermalTrees/a4all_%d-%d_max_30001_sinsub_10_3_ad_2.root",i,i+39);
      TString adding = TString::Format("thermalTrees/simulated_%d-%d_max_1000_sinsub_10_3_ad_2.root",i,i+39);
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
      sumfile = new TFile(TString::Format("/Volumes/SDCard/data/%s/%d_max_30001_sinsub_10_3_ad_2.root", "a4all",run));; 
      gROOT->cd(); 

      if (eventsfile) delete eventsfile; 
      eventsfile = new TFile(TString::Format("source_maps/%d_%d.root", run,run)); 
      gROOT->cd(); 

      sumtree = (TTree*) sumfile->Get("anita4"); 
      eventstree = (TTree*) eventsfile->Get("events"); 

      sumtree->SetBranchAddress("summary",&sum); 
      sumtree->SetBranchAddress("pat",&gps); 
      sumtree->BuildIndex("eventNumber"); 
      eventstree->BuildIndex("event"); 
      eventstree->SetBranchAddress("S",&S); 

      loaded_run =run; 
    }

    int event = events[i]; 
    int iteration = iters[i]; 
    sumtree->GetEntryWithIndex(event); 
    eventstree->GetEntryWithIndex(event); 
    printf("%d %d %d %g\n", run,event, iteration, S); 
    m->add(sum,gps,AnitaPol::AnitaPol_t (iteration / 5), iteration % 5, -S); 
  }
  
  f.cd(); 
  m->Write(out_key); 
  return events.size(); 
}





int evaluateSourceMap(int start_run = 50, int end_run = 367, 
                      bool mc = false, bool decimated = true,
                      const char * prefix = "source_maps_eval/",
                      const char * infile="all_source_maps.root",
                      const char * key = "map_weighted_culled", 
                      const char * removed_events_file = "removed_events.txt", 
                      bool weighted = true, bool vpol_only = false ) 
{
  TChain c(mc ? "simulation":"anita4"); 
  if (mc) 
  {
//    c.Add("thermalTrees/simulated*.root"); 
    c.Add("thermalTrees/simulated*.root"); 
    decimated = false; 
  }
  else
  {
    addRuns(c,start_run,end_run); 
  }

  TFile of(TString::Format("%s%d_%d_%s.root",prefix,start_run,end_run, decimated ? "10pct" : mc ? "sim" : "full" ),"RECREATE"); 

  TTree * ot = new TTree("overlap","Overlap"); 
  double O; 
  double S; 
  int run;
  int event;
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


  ot->Branch("O",&O); 
  ot->Branch("S",&S); 
  ot->Branch("run",&run); 
  ot->Branch("F",&F); 
  ot->Branch("event",&event); 
  ot->Branch("pol",&pol); 
  ot->Branch("theta",&theta); 
  ot->Branch("nclustered",&nclustered,"nclustered[10]/I"); 
  ot->Branch("polangle",&polangle); 
  ot->Branch("base_sum",&base_sum); 
  ot->Branch("max_base_index",&max_base_index); 
  ot->Branch("max_base_p",&max_base_p); 
  ot->Branch("weight",&wgt); 
  ot->Branch("mcE",&mcE); 
  ot->Branch("peak",&peak); 
  ot->Branch("removed",&removed); 

  TFile f(infile); 
  UCorrelator::ProbabilityMap * map = (UCorrelator::ProbabilityMap*) f.Get(key); 
  UCorrelator::ProbabilityMap * source_map = map; //for mc 
  UCorrelator::ProbabilityMap::Params * map_pars = map_params(); 
  int n = c.Draw("run:eventNumber:iteration:F",  TCut(TString::Format("(%s) && (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 

  std::vector<std::vector<double> > counts(10, std::vector<double> (map->segmentationScheme()->NSegments())); 

  for (int d = 0; d < 10; d++)
  {
    map->groupAdjacent(map->getWgtAboveLevel(d), 0, &counts[d][0]); 
  }


  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  std::set<int> * removed_events = removed_events_file ? getRemovedEvents(removed_events_file) : 0; 

  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  AnitaDataset * d = 0; 
  for (int i = 0; i < n; i++) 
  {
    run = c.GetV1()[i]; 
    if (run!= loaded_run)
    {
        if (sumfile) delete sumfile; 
        if (d) delete d; 

        if (decimated) 
        {
          d = new AnitaDataset(run,true); 
        }

        sumfile = new TFile(TString::Format("/Volumes/SDCard/data/%s/%d_max_30001_sinsub_10_3_ad_2.root", mc ? "simulated" : "a4all",run));; 
        gROOT->cd(); 
        sumtree = (TTree*) sumfile->Get(mc ? "simulation" : "anita4"); 
        sumtree->SetBranchAddress("summary",&sum); 
        sumtree->SetBranchAddress("pat",&gps); 
        sumtree->BuildIndex("eventNumber"); 
        loaded_run =run; 
    }
      
    event = int(c.GetV2()[i]); 


    if (decimated)
    {

      int ret = d->getEvent(event,true); 
      printf("d->getEvent(%d) returns: %d\n", event, ret); 
      if (ret < 0) continue; 
    }

    S = c.GetW()[i]; 
    sumtree->GetEntryWithIndex(event); 
    pol = int(c.GetV3()[i]) / 5; 
    peak = int(c.GetV3()[i]) % 5; 
    polangle = 90/TMath::Pi() * TMath::ATan2(sum->deconvolved_filtered[pol][peak].max_dU, sum->deconvolved_filtered[pol][peak].max_dQ); 
    F = c.GetV4()[i]; 
    theta = sum->peak[pol][peak].theta; 

    removed = (removed_events && removed_events->count(event));  

    double phi = sum->peak[pol][peak].phi; 

    if (mc) 
    {
      wgt = sum->mc.weight; 
      mcE = sum->mc.energy; 
      if ( fabs(sum->mc.theta-theta) > 4) continue; 
      if ( fabs(FFTtools::wrap(sum->mc.phi-phi, 360,0)) > 4) continue; 
    }


    if ( (!weighted && F < cutoff) || (vpol_only && pol == 0)) continue; 

    std::vector<std::pair<int,double> > bases;
    std::vector<std::pair<int,double> > dens;
    double inv_two_pi_sqrt_det; 

    if (mc) 
    {
      map = new UCorrelator::ProbabilityMap(map_pars); 
      map->add(sum,gps,AnitaPol::AnitaPol_t(pol),peak,S); 
      map->combineWith(*source_map); 

      for (int d = 0; d< 10; d++) 
      {
        
        map->groupAdjacent(map->getWgtAboveLevel(d), 0, &counts[d][0]); 
      }

    }

    O = map->overlap(sum,gps,AnitaPol::AnitaPol_t(pol),peak,true,S, &bases, UCorrelator::ProbabilityMap::OVERLAP_SUM_SQRTS, !removed,0,&dens,&inv_two_pi_sqrt_det) / sqrt(S); 

    for (int d = 0; d < 10; d++) 
    {
      nclustered[d] = -1; 
      double thresh = UCorrelator::ProbabilityMap::dist2dens( map->getLevel(d), inv_two_pi_sqrt_det); 

      for (unsigned si = 0; si < dens.size(); si++)
      {
        if (dens[si].second < thresh) continue; 

        if (counts[d][dens[si].first])
        {
           nclustered[d] = (int) counts[d][dens[si].first]; 
           break; 
        }
      }
    }

    if (mc) delete map; 
   
//    if (O < 0) O = -1; 
    printf("r/e: %d/%d/%g/%g\n",run,event, F, O); 

    max_base_index = -1;
    base_sum = 0;
    max_base_p = 0;
    for (unsigned i = 0; i < bases.size();i++)
    {
      if (bases[i].second > max_base_p)
      {
        max_base_index = bases[i].first;
        max_base_p = bases[i].second; 
      }
      base_sum += bases[i].second; 
    }

    of.cd(); 
    ot->Fill(); 
  }

  of.cd(); 
  ot->Write(); 
  return 0; 

}





// int makeSourceMap(int start_run = 50, int end_run = 367, const char * prefix = "source_maps/test_")
int makeSourceMap(int start_run = 1, int end_run = 200, const char * prefix = "source_maps/test_")
{

  // Start getting the run / event numbers of events that pass our cut

  // TChain c("anita4"); 
  TChain c("simulation"); 

  addRuns(c,start_run,end_run); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  int n = c.Draw("run:eventNumber:iteration:F",  TCut(TString::Format("(%s) && (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 
  printf("%d events pass selection\n", n); 


  UCorrelator::ProbabilityMap::Params * p = map_params(); 

  TFile * f  = new TFile(TString::Format("%s%d_%d.root",prefix,start_run, end_run), "RECREATE"); 
  UCorrelator::ProbabilityMap *map_weighted = new UCorrelator::ProbabilityMap(p); 
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
  for (int i = 0; i < 1; i+=1) 
  {
    run = c.GetV1()[i]; 
    if (run!= loaded_run)
    {
        if (sumfile) delete sumfile; 
        // sumfile = new TFile(TString::Format("/Volumes/SDCard/data/%s/%d_max_30001_sinsub_10_3_ad_2.root", "a4all",run));; 
        sumfile = new TFile(TString::Format("/Volumes/SDCard/data/%s/%d_max_1000_sinsub_10_3_ad_2.root", "simulated",run));; 
        gROOT->cd(); 
        // sumtree = (TTree*) sumfile->Get("anita4"); 
        sumtree = (TTree*) sumfile->Get("simulation"); 
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

    nsegs = map_weighted->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, S); 

    double integ = map_weighted->getProbSumsIntegral(false); 
    double integ_norm = map_weighted->getProbSumsIntegral(true); 
    dinteg = integ-last_integ; 
    dinteg_norm = integ_norm-last_integ_norm; 
    last_integ_norm = integ_norm; 
    last_integ = integ_norm; 

    tr->Fill(); 

    if (F > cutoff) 
    {
      map->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, 1); 
    }


  }


  f->cd(); 
  map->Write("map_unweighted"); 
  map_weighted->Write("map_weighted"); 
  tr->Write(); 
  return 0; 
}




