#include "FFTtools.h" 
#include "AnitaConventions.h" 
#include "AnitaDataset.h" 

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



double cutoff = 3; 
const char * weight = "((F > 3.25) + (F < 3.25 && F > 2) * exp (-((abs(F-3.25))^0.5879) / 0.4231 )) * ( F > 0 && theta > 3 && isMostImpulsive && !payloadBlast && MaxPeak < 1000 && chisq < 0.01 && theta < 40)";


void addRuns(TChain & c, int start_run, int end_run)
{
  for (int i = 130; i < 430; i+=10) 
  {
    if (start_run < i + 10 && end_run >= i)
    {
      c.Add(TString::Format("thermalTrees/a3all_%d-%d_sinsub_10_3_ad_2.root",i,i+9)); 
    }
  }
}


int evaluateSourceMap(int start_run = 300, int end_run = 330, 
                      bool mc = false, bool decimated = true,
                      const char * prefix = "source_maps_eval/",
                      const char * infile="all_source_maps.root",
                      const char * key = "map_weighted", 
                      bool weighted = false, bool vpol_only = true ) 
{
  TChain c(mc ? "simulation":"anita3"); 
  if (mc) 
  {
    c.Add("thermalTrees/simulated*.root"); 
    decimated = false; 
  }
  else
  {
    //add the necessary trees 

    for (int i = 130; i < 430; i+=10) 
    {
      if (start_run < i + 10 && end_run >= i)
      {
        c.Add(TString::Format("thermalTrees/a3all_%d-%d_sinsub_10_3_ad_2.root",i,i+9)); 
      }
    }
  }

  TFile of(TString::Format("%s%d_%d_%s.root",prefix,start_run,end_run, decimated ? "10pct" : mc ? "sim" : "full" ),"RECREATE"); 

  TTree * ot = new TTree("overlap","Overlap"); 
  double O; 
  double S; 
  int run;
  int event;
  double theta; 
  double base_sum; 
  double F; 
  double wgt = 1; 
  int max_base; 
  
  ot->Branch("O",&O); 
  ot->Branch("S",&S); 
  ot->Branch("run",&run); 
  ot->Branch("F",&F); 
  ot->Branch("event",&event); 
  ot->Branch("theta",&theta); 
  ot->Branch("base_sum",&base_sum); 
  ot->Branch("max_base",&max_base); 
  ot->Branch("weight",&wgt); 

  TFile f(infile); 
  UCorrelator::ProbabilityMap * map = (UCorrelator::ProbabilityMap*) f.Get(key); 
  int n = c.Draw("run:eventNumber:iteration:F",  TCut(TString::Format("(%s) * (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

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

        sumfile = new TFile(TString::Format("%s/%d_sinsub_10_3_ad_2.root", mc ? "simulated" : "a3all",run));; 
        gROOT->cd(); 
        sumtree = (TTree*) sumfile->Get(mc ? "simulation" : "anita3"); 

        sumtree->SetBranchAddress("summary",&sum); 
        sumtree->SetBranchAddress("pat",&gps); 
        sumtree->BuildIndex("eventNumber"); 
        loaded_run =run; 
    }
      
    event = int(c.GetV2()[i]); 

    if (decimated && d->getEvent(event,true) < 0) continue; 

    S = c.GetW()[i]; 
    sumtree->GetEntryWithIndex(event); 
    int pol = int(c.GetV3()[i]) / 5; 
    int peak = int(c.GetV3()[i]) % 5; 
    F = c.GetV4()[i]; 
    theta = sum->peak[pol][peak].theta; 
    double phi = sum->peak[pol][peak].phi; 

    if (mc) 
    {
      wgt = sum->mc.weight; 
      if ( fabs(sum->mc.theta-theta) > 4) continue; 
      if ( fabs(FFTtools::wrap(sum->mc.phi-phi, 360,0)) > 4) continue; 
    }


    if ( (!weighted && F < cutoff) || (vpol_only && pol == 0)) continue; 

    std::vector<std::pair<int,double> > bases;
    O = map->overlap(sum,gps,AnitaPol::AnitaPol_t(pol),peak,true,S, &bases, UCorrelator::ProbabilityMap::OVERLAP_SUM_SQRTS, !mc) / S; 
    printf("r/e: %d/%d/%g/%g\n",run,event, F, O); 

    max_base = -1;
    base_sum = 0;
    double base_max_p = 0;
    for (unsigned i = 0; i < bases.size();i++)
    {
      if (bases[i].second > base_max_p)
      {
        max_base = bases[i].first;
        base_max_p = bases[i].second; 
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





int makeSourceMap(int start_run = 300, int end_run = 330, const char * prefix = "source_maps/")
{

  // Start getting the run / event numbers of events that pass our cut

  TChain c("anita3"); 

  addRuns(c,start_run,end_run); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  int n = c.Draw("run:eventNumber:iteration:F",  TCut(TString::Format("(%s) * (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 

  StereographicGrid g(1024,1024); 

  TF1 f_dtheta("ftheta", "[0] / x^[1]", 1, 50);
  TF1 f_dphi("fphi", "[0] / x^[1]", 1, 50);
  f_dtheta.SetParameter(0, 0.3936); 
  f_dtheta.SetParameter(1, 0.2102); 
  f_dphi.SetParameter(0, 1.065); 
  f_dphi.SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel m1(f_dtheta, f_dphi, true,true);
  UCorrelator::PointingResolutionModelPlusHeadingError m(20, &m1); 

  Refraction::SphRay ref; 

  UCorrelator::ProbabilityMap::Params p; 
  p.refract = &ref; 
  p.seg = &g; 
  p.point = &m; 
  p.collision_detection = false; 
 

  TFile * f  = new TFile(TString::Format("%s%d_%d.root",prefix,start_run, end_run), "RECREATE"); 
  UCorrelator::ProbabilityMap *map_weighted = new UCorrelator::ProbabilityMap(&p); 
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(&p); 

  int loaded_run = 0; 
  TFile * sumfile = 0; 
  TTree * sumtree = 0; 
  for (int i = 0; i < n; i++) 
  {
    int run = c.GetV1()[i]; 
    if (run!= loaded_run)
    {
        if (sumfile) delete sumfile; 
        sumfile = new TFile(TString::Format("%s/%d_sinsub_10_3_ad_2.root", "a3all",run));; 
        gROOT->cd(); 
        sumtree = (TTree*) sumfile->Get("anita3"); 
        sumtree->SetBranchAddress("summary",&sum); 
        sumtree->SetBranchAddress("pat",&gps); 
        sumtree->BuildIndex("eventNumber"); 
        loaded_run = run; 
    }
      
    int ev = int(c.GetV2()[i]); 
    double S = c.GetW()[i]; 
    sumtree->GetEntryWithIndex(ev); 
    int pol = int(c.GetV3()[i]) / 5; 
    int peak = int(c.GetV3()[i]) % 5; 
    double F = c.GetV4()[i]; 
    printf("r/e: %d/%d/%g(%g)\n",run,ev, F, S); 
    map_weighted->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, S); 
    if (F > cutoff) 
    {
      map->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, 1); 
    }
  }


  f->cd(); 
  map->Write("map_unweighted"); 
  map_weighted->Write("map_weighted"); 
  return 0; 
}
