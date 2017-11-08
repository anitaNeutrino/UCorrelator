int makeHotspotMap(int run, int max = 100, int start = 0,  const char * fmt = "a3all/%d_sinsub_10_3_ad_2.root", const char * tree = "anita3") 
{
  TChain c(tree); 
  c.Add(TString::Format(fmt, run)); 

  Refraction::SphRay ref; 
  StereographicGrid *g = new StereographicGrid(1024,1024); 
   //my formula
  TF1 f_dtheta("ftheta", "[0] / x^[1]", 1, 50);
  TF1 f_dphi("fphi", "[0] / x^[1]", 1, 50);
  f_dtheta.SetParameter(0, 0.3936); 
  f_dtheta.SetParameter(1, 0.2102); 
  f_dphi.SetParameter(0, 1.065); 
  f_dphi.SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel m(f_dtheta, f_dphi, true);

  UCorrelator::ProbabilityMap::Params p; 
  p.seg = g; 
  p.refract = &ref; 
  p.point = &m; 
  p.collision_detection = false; 
  p.backwards_params.el_cutoff = 0; 
  p.backwards_params.num_samples_per_bin = 9; 
  p.maximum_distance=20; 


  AnitaEventSummary * sum = 0;; 
  Adu5Pat * pat = 0;; 

  c.SetBranchAddress("summary",&sum); 
  c.SetBranchAddress("pat",&pat); 

  TFile f(TString::Format("hotspot/%d_%d_%d.root",run,max,start),"RECREATE"); 

  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(&p); 

  for (int i = start; i < (max > 0 ? start+max : c.GetEntries()) ; i++) 
  {
    c.GetEntry(i); 
    printf("%d\n",i); 
    for (int pol = 0; pol < 2; pol++)
    {
      for (int peak = 0; peak < 1; peak++) 
      {
        if (sum->peak[pol][peak].theta > 5)
        {
          map->add(sum, pat, AnitaPol::AnitaPol_t(pol), peak, 1, max == 1 ? &f: 0 ); 
        }
      }
    }
  }

  map->segmentationScheme()->Draw("colz",map->getProbSums(true)); 
  f.cd(); 
  map->Write("hotspot"); 

  return 0; 
}
