
#include "FFTtools.h" 
#include "AnitaConventions.h" 

UCorrelator::ProbabilityMap* testMapWais(int run =342, int max = 20, int nskip = 0) 
{

  FFTtools::loadWisdom("wisdom.dat"); 
  AnitaDataset d(run); 


  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 

  UCorrelator::Analyzer  analyzer (&cfg); 
  AnitaEventSummary sum; 
  FilterStrategy strategy; 
  UCorrelator::fillStrategyWithKey(&strategy,"sinsub_10_3_ad_2"); 

  TF1 f_dtheta("ftheta", "[0] / x^[1]", 1, 50);
  TF1 f_dphi("fphi", "[0] / x^[1]", 1, 50);
  f_dtheta.SetParameter(0, 0.3936); 
  f_dtheta.SetParameter(1, 0.2102); 

  f_dphi.SetParameter(0, 1.065); 
  f_dphi.SetParameter(1, 0.2479); 

  UCorrelator::PointingResolutionParSNRModel m(&f_dtheta, &f_dphi, true);
 
  StereographicGrid *g = new StereographicGrid(4096,4096); 
  UCorrelator::ProbabilityMap::Params p; 
  Refraction::SphRay ref; 
  p.refract = &ref; 
  p.seg = g; 
  p.point = &m; 
  p.collision_detection = true; 
  p.mc_params.n = 1e5; 
  p.projection = UCorrelator::ProbabilityMap::Params::MC; 
  p.backwards_params.el_cutoff = 0; 

  UCorrelator::ProbabilityMap::Params p2; 
  p2.seg = g; 
  p2.refract = &ref; 
  p2.point = &m; 
  p2.collision_detection = true; 
  p2.backwards_params.el_cutoff = 0; 
 
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(&p); 
  UCorrelator::ProbabilityMap *map2 = new UCorrelator::ProbabilityMap(&p2); 
  int ndone = 0; 

  TGraph *  g_anita = new TGraph; 
  TGraph *  g_proj = new TGraph; 


  for (int i = 0; i< d.N(); i++)
  {
    d.getEntry(i); 
    printf("----(%d)-----\n",i); 
    UsefulAdu5Pat pat(d.gps()); 
    if (UCorrelator::isWAISHPol(&pat, d.header()))
    {
      if (nskip-- > 0) continue; 
      printf("Processing event %d (%d)\n",d.header()->eventNumber,ndone); 
      FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 
      analyzer.analyze(&ev, &sum); 
//      printf("dphi: %g; dtheta: %g, theta: %g, theta_adjust: %g\n", FFTtools::wrap(sum.peak[0][0].phi - sum.wais.phi, 360,0), FFTtools::wrap(sum.peak[0][0].theta - sum.wais.theta,360,0), sum.peak[0][0].theta, sum.peak[0][0].theta_adjustment_needed); 
      AntarcticCoord a(AntarcticCoord::WGS84, d.gps()->latitude, d.gps()->longitude, d.gps()->altitude); 
      a.to(AntarcticCoord::STEREOGRAPHIC); 
      g_anita->SetPoint(g_anita->GetN(),  a.x, a.y); 

      PayloadParameters project; 
      int status = PayloadParameters::findSourceOnContinent(sum.peak[0][0].theta, sum.peak[0][0].phi, d.gps(), & project, &ref, p.collision_params.dx); 
      project.source.to(AntarcticCoord::STEREOGRAPHIC); 
      printf("%d %g %g %g\n", status, project.source.x, project.source.y, project.source.z); 
      g_proj->SetPoint(g_proj->GetN(),  project.source.x, project.source.y); 


  //    map->add(&sum,d.gps(), AnitaPol::kHorizontal, 0); 
      map2->add(&sum,d.gps(), AnitaPol::kHorizontal, 0); 
      map->dumpNonZeroBases(); 

      ndone++;
    }


    if (max && ndone >= max) break; 

  }

  FFTtools::saveWisdom("wisdom.dat"); 
  
  map->dumpNonZeroBases(); 

  //draw antarctica 
//  TProfile2D * bg = RampdemReader::getMap(RampdemReader::surface, 10); 
//  bg->SetMaximum(1000); 
//  bg->DrawCopy("cont"); 

  //draw the p map 
  gStyle->SetPalette(kTemperatureMap); 

  //draw WAIS 
  double wais_lon=AnitaLocations::getWaisLongitude(); 
  double wais_lat=AnitaLocations::getWaisLatitude(); 
  double wais_E, wais_N; 
  RampdemReader::LonLatToEastingNorthing(wais_lon,wais_lat, wais_E, wais_N); 
  double range[4]; 
  range[0] = wais_E - 250e3; 
  range[1] = wais_E + 300e3; 
  range[2] = wais_N - 200e3; 
  range[3] = wais_N + 200e3; 

  TMarker * wais = new TMarker(wais_E, wais_N, 3); 
  g_anita->SetLineColor(2); 
  g_anita->SetMarkerColor(2); 
  g_anita->SetMarkerStyle(4); 

  g_proj->SetLineColor(3); 
  g_proj->SetMarkerColor(6); 
  g_proj->SetMarkerStyle(5); 

  TCanvas *c = new TCanvas; 
  c->Divide(2,2); 
  c->cd(1)->SetLogz(); 
  map->segmentationScheme()->Draw("colz", map2->getDensitySums(),range); 
  wais->SetMarkerColor(3); 
  wais->Draw("psame"); 
  g_anita->Draw("lpsame"); 
  g_proj->Draw("psame"); 

  c->cd(2)->SetLogz(); 
  map->segmentationScheme()->Draw("colz", map2->getDensitySumsNormalized(),range); 
  g_anita->Draw("lpsame"); 
  g_proj->Draw("psame"); 
  wais->Draw("psame"); 

  c->cd(3); 

  map->segmentationScheme()->DrawI("colz", map2->getNAboveLevel(0),range); 
  g_proj->Draw("psame"); 
  wais->Draw("psame"); 
  /*
  std::vector<double> ratio (g.NSegments()); 

  for (size_t i = 0; i < ratio.size(); i++)
  {

    if (map->getDensitySums()[i] && map2->getDensitySums()[i]) 
    {
      ratio[i] = (map->getDensitySums()[i]  - map2->getDensitySums()[i]) / map2->getDensitySums()[i]; 
    }

  }
  */
//  map->segmentationScheme()->Draw("colz", &ratio[0],range); 
  g_anita->Draw("lpsame"); 
  wais->Draw("psame"); 
  c->cd(4); 

  TH2 * huge = new TH2D("surface","surface", 1024, range[0], range[1], 1024,range[2],range[3]); 

  for (int i = 1; i<= huge->GetNbinsX(); i++)
  {
    for (int j = 1; j <= huge->GetNbinsY(); j++)
    {
      huge->SetBinContent(i,j, RampdemReader::SurfaceAboveGeoidEN( huge->GetXaxis()->GetBinCenter(i), huge->GetYaxis()->GetBinCenter(j), RampdemReader::rampdem)); 
    }
  }

  huge->SetStats(0); 
  huge->Draw("colz"); 
  wais->Draw("psame"); 

  map->dumpNonZeroBases(); 

  printf("wais: %g\n", map->getBaseDensitySums()[293]); 

  map->dumpNonZeroBases(); 

  return map; 
}
