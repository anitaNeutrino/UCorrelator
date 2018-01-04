
#include "AnitaConventions.h" 
#include "AnitaEventSummary.h" 

void testMapWaisOverlap(int run0 =342, int ev0 = 0, int run1=343, int ev1 = 0, bool collision_detect = true,  const char * filter ="sinsub_10_3_ad_2") 
{

  TChain c0("wais"); 
  TChain c1("wais"); 

  c0.Add(TString::Format("wais/%d_%s.root", run0, filter)); 
  c1.Add(TString::Format("wais/%d_%s.root", run1, filter)); 

  //my formula
//  TF1 f_dtheta("ftheta", "[0] / x^[1]", 1, 50);
//  TF1 f_dphi("fphi", "[0] / x^[1]", 1, 50);
//  f_dtheta.SetParameter(0, 0.3936); 
//  f_dtheta.SetParameter(1, 0.2102); 
//  f_dphi.SetParameter(0, 1.065); 
//  f_dphi.SetParameter(1, 0.2479); 
//  UCorrelator::PointingResolutionParSNRModel m(&f_dtheta, &f_dphi, true);

  // Ben R's formula
  TF1 sigmaPhi("sigmaPhi","[0]/(x) + [1]",4,100);
  sigmaPhi.SetParameter(0,0.68726);
  sigmaPhi.SetParameter(1,0.338725);
       
  TF1 sigmaTheta("sigmaTheta","[0]/(x) + [1]",4,100);
  sigmaTheta.SetParameter(0,0.132516);
  sigmaTheta.SetParameter(1,0.1557);
  UCorrelator::PointingResolutionParSNRModel m(sigmaTheta, sigmaPhi, false);



  Refraction::SphRay ref; 
  StereographicGrid *g = new StereographicGrid(1024,1024); 
  UCorrelator::ProbabilityMap::Params p; 
  p.seg = g; 
  p.refract = &ref; 
  p.point = &m; 
  p.collision_detection = collision_detect; 
  p.backwards_params.el_cutoff = 0; 
  p.backwards_params.num_samples_per_bin = 16; 
  p.maximum_distance=20; 


  AnitaEventSummary * sum0 = 0;; 
  AnitaEventSummary * sum1 = 0;; 

  Adu5Pat * gps0 = 0; 
  Adu5Pat * gps1 = 0; 

  c0.SetBranchAddress("summary",&sum0); 
  c1.SetBranchAddress("summary",&sum1); 
  c0.SetBranchAddress("pat",&gps0); 
  c1.SetBranchAddress("pat",&gps1); 

  c0.BuildIndex("eventNumber"); 
  c1.BuildIndex("eventNumber"); 
  ev0 >0 ? c0.GetEntryWithIndex(ev0) : c0.GetEntry(-ev0); 
  ev1 >0 ? c1.GetEntryWithIndex(ev1) : c1.GetEntry(-ev1); 
  
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(&p); 

  TFile f("overlap.root","RECREATE"); 
  f.mkdir("triangles"); 
  map->add(sum0, gps0, AnitaPol::kHorizontal,0, 1,&f); 
  map->add(sum1, gps1, AnitaPol::kHorizontal,0); 


  AntarcticCoord a0(AntarcticCoord::WGS84, gps0->latitude, gps0->longitude, gps0->altitude); 
  AntarcticCoord a1(AntarcticCoord::WGS84, gps1->latitude, gps1->longitude, gps1->altitude); 

  AntarcticCoord p0(AntarcticCoord::WGS84, sum0->peak[0][0].latitude, sum0->peak[0][0].longitude, sum0->peak[0][0].altitude); 
  AntarcticCoord p1(AntarcticCoord::WGS84, sum1->peak[0][0].latitude, sum1->peak[0][0].longitude, sum1->peak[0][0].altitude); 

  printf(" (%g,%g,%g), (%g,%g,%g) \n", p0.x,p0.y,p0.z, p1.x,p1.y,p1.z); 

  double wais_lon=AnitaLocations::getWaisLongitude(); 
  double wais_lat=AnitaLocations::getWaisLatitude(); 
  double wais_E, wais_N; 
  RampdemReader::LonLatToEastingNorthing(wais_lon,wais_lat, wais_E, wais_N); 
  TMarker * wais = new TMarker(wais_E, wais_N, 3); 
  wais->SetMarkerColor(4); 
 

  p0.to(AntarcticCoord::STEREOGRAPHIC); 
  p1.to(AntarcticCoord::STEREOGRAPHIC); 
  a0.to(AntarcticCoord::STEREOGRAPHIC); 
  a1.to(AntarcticCoord::STEREOGRAPHIC); 

  TMarker * m0 = new TMarker(p0.x, p0.y, 3); 
  TMarker * m1 = new TMarker(p1.x, p1.y, 3); 
  TArrow * ar0 = new TArrow(a0.x,a0.y,p0.x,p0.y); 
  TArrow * ar1 = new TArrow(a1.x,a1.y,p1.x,p1.y); 

  m0->SetMarkerColor(3); 
  m1->SetMarkerColor(6); 
  ar0->SetLineColor(3); 
  ar1->SetLineColor(6); 

  map->segmentationScheme()->Draw("colz", map->getProbSums(true)); 

  m0->Draw("psame"); 
  m1->Draw("psame"); 
  ar0->Draw("psame"); 
  ar1->Draw("psame"); 
  wais->Draw("psame"); 

  double overlap0 = map->overlap(sum0,gps0,AnitaPol::kHorizontal,0,false,1,0, UCorrelator::ProbabilityMap::OVERLAP_MAXES); 
  double overlap1 = map->overlap(sum1,gps1,AnitaPol::kHorizontal,0,false,1,0); 

  printf("%g %g\n", overlap0, overlap1); 

  overlap0 = map->overlap(sum0,gps0,AnitaPol::kHorizontal,0,true,1,0); 
  overlap1 = map->overlap(sum1,gps1,AnitaPol::kHorizontal,0,true,1,0); 

  printf("%g %g\n", overlap0, overlap1); 


} 
