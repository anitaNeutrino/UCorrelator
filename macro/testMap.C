#include "AnitaConventions.h" 
#include "AnitaEventSummary.h" 
#include "AnitaDataset.h" 

UCorrelator::ProbabilityMap* testMap(int run, int event_number, int pol = 1, int peak = 0) 
{

  TChain c("anita3"); 
  c.Add(TString::Format("a3all/%d*.root",run)); 
  c.BuildIndex("eventNumber"); 
  AnitaEventSummary * sum = 0; 
  Adu5Pat * pat = 0; 
  c.SetBranchAddress("summary",&sum); 
  c.SetBranchAddress("pat",&pat); 

  c.GetEntryWithIndex(event_number); 

  TFile f("debug_data.root","RECREATE"); 

  StereographicGrid  * g = new StereographicGrid(1024,1024); 

  TF1 f_dtheta("ftheta", "[0] / x^[1]", 1, 50);
  TF1 f_dphi("fphi", "[0] / x^[1]", 1, 50);
  f_dtheta.SetParameter(0, 0.3936); 
  f_dtheta.SetParameter(1, 0.2102); 
  f_dphi.SetParameter(0, 1.065); 
  f_dphi.SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel * m1 = new UCorrelator::PointingResolutionParSNRModel(f_dtheta, f_dphi, true,true);
  UCorrelator::PointingResolutionModelPlusHeadingError*  m = new UCorrelator::PointingResolutionModelPlusHeadingError(20, m1); 

  Refraction::SphRay * ref = new Refraction::SphRay; 

  UCorrelator::ProbabilityMap::Params p; 
  p.refract = ref; 
  p.seg = g; 
  p.point = m; 
  p.collision_detection = false; 
  p.backwards_params.enhance_threshold = 1e-2; 
//  p.backwards_params.num_samples_per_bin=16;
  p.verbosity = 3;

  printf("THETA: %g\n", sum->peak[pol][peak].theta); 
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(&p); 
  map->add( sum, pat, AnitaPol::AnitaPol_t(pol), peak,1,&f); 
  map->segmentationScheme()->Draw("colz", map->getProbSums(true)); 

  f.cd(); 
  map->Write("map"); 

  return map ; 

}
