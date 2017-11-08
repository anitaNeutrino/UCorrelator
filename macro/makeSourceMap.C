#include "FFTtools.h" 
#include "AnitaConventions.h" 

int makeSourceMap(int start_run = 300, int end_run = 330 , const char * weight = "( (Fisher > 5) + (Fisher < 5 ) * 0.8864*TMath::Gaus(Fisher,5,0.95)) *  ( Fisher > 0 && theta > 0 && isMostImpulsive  && !isTraining && !payloadBlast && !payloadBlastExtended && blastFraction < 0.2 && nSectorsWhereBottomExceedsTop < 28 && nSectorsWhereBottomExceedsTop > 4)", bool mc= false) 
{

  // Start getting the run / event numbers of events that pass our cut

  TChain c(mc ? "simulation":"anita3"); 
  if (mc) 
  {
    c.Add("thermalTrees/simulated*.root"); 
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


  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * gps = new Adu5Pat; 

  int n = c.Draw("run:eventNumber:iteration",  TCut(TString::Format("(%s) * (run >= %d && run <= %d)", weight,start_run,end_run)) , "goff" ); 

  StereographicGrid g(1024,1024); 

  TF1 f_dtheta("ftheta", "[0] / x^[1]", 1, 50);
  TF1 f_dphi("fphi", "[0] / x^[1]", 1, 50);
  f_dtheta.SetParameter(0, 0.3936); 
  f_dtheta.SetParameter(1, 0.2102); 
  f_dphi.SetParameter(0, 1.065); 
  f_dphi.SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel m1(f_dtheta, f_dphi, true,true);
  UCorrelator::PointingResolutionModelPlusHeadingError m(60, &m1); 



  Refraction::SphRay ref; 

  UCorrelator::ProbabilityMap::Params p; 
  p.refract = &ref; 
  p.seg = &g; 
  p.point = &m; 
  p.collision_detection = false; 
 
  TFile * f  = new TFile(TString::Format("source_map_%d_%d.root",start_run, end_run), "RECREATE"); 
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
        sumfile = new TFile(TString::Format("%s/%d_sinsub_10_3_ad_2.root", mc ? "simulated" : "a3all",run));; 
        gROOT->cd(); 
        sumtree = (TTree*) sumfile->Get(mc ? "simulation" : "anita3"); 

        sumtree->SetBranchAddress("summary",&sum); 
        sumtree->SetBranchAddress("pat",&gps); 
        sumtree->BuildIndex("eventNumber"); 
    }
      
    int ev = int(c.GetV2()[i]); 
    printf("r/e: %d/%d\n",run,ev); 
    sumtree->GetEntryWithIndex(ev); 
    int pol = int(c.GetV3()[i]) / 5; 
    int peak = int(c.GetV3()[i]) % 5; 
    map->add(sum, gps, AnitaPol::AnitaPol_t(pol), peak, c.GetW()[i]); 
  }


  f->cd(); 
  map->Write(); 
  return 0; 
}
