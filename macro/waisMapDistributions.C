
double get_sum (int N, const double *v)
{
  double ans = 0; 
  for (int i = 0; i < N; i++) ans+= v[i]; 
  return ans; 
}


int waisMapDistributions(int N = 1000, int iout = 0, const char * output = "waismap/WaisMap_%d.root")
{

  TChain c("wais"); 
  c.Add("wais/*.root"); 

  int Nwais = c.GetEntries(); 
  AnitaEventSummary * sum = 0; 
  Adu5Pat * pat = 0; 
  c.SetBranchAddress("summary",&sum); 
  c.SetBranchAddress("pat",&pat); 



  gRandom->SetSeed(iout); 


  //my formula
  TF1 f_dtheta("ftheta", "[0] / x^[1]", 1, 50);
  TF1 f_dphi("fphi", "[0] / x^[1]", 1, 50);
  f_dtheta.SetParameter(0, 0.3936); 
  f_dtheta.SetParameter(1, 0.2102); 
  f_dphi.SetParameter(0, 1.065); 
  f_dphi.SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel m1(f_dtheta, f_dphi, true,true);
  UCorrelator::PointingResolutionModelPlusHeadingError m(20, &m1); 



  // Ben R's formula
//  TF1 sigmaPhi("sigmaPhi","[0]/(x) + [1]",4,100);
//  sigmaPhi.SetParameter(0,0.68726);
//  sigmaPhi.SetParameter(1,0.338725);
//       
//  TF1 sigmaTheta("sigmaTheta","[0]/(x) + [1]",4,100);
//  sigmaTheta.SetParameter(0,0.132516);
//  sigmaTheta.SetParameter(1,0.1557);
//  UCorrelator::PointingResolutionParSNRModel m(&sigmaTheta, &sigmaPhi, false);

  Refraction::SphRay ref; 
  StereographicGrid *g = new StereographicGrid(1024,1024); 
  UCorrelator::ProbabilityMap::Params p; 
  p.seg = g; 
  p.refract = &ref; 
  p.point = &m; 
  p.collision_detection = false; 
  p.backwards_params.el_cutoff = 0; 

 
  TFile fout(TString::Format(output,iout),"RECREATE"); 
  TTree * tree = new TTree("waisMap","WAIS Map Distributions"); 

  double overlap_12; 
  double overlap_12_norm; 
  double wais_density_1; 
  double wais_density_2; 
  double wais_distance_1; 
  double wais_distance_2; 
  double wais_density_1_norm; 
  double wais_density_2_norm; 
  double wais_distance_1_norm; 
  double wais_distance_2_norm; 
  double sum_1, sum_2; 
  double dphi1, dphi2; 
  double dtheta1, dtheta2; 
  double snr_1;
  double snr_2;
  double L; 
  int run_1; 
  int run_2; 
  int event_1; 
  int event_2; 

  tree->Branch("run_1",&run_1);
  tree->Branch("run_2",&run_2);
  tree->Branch("event_1",&event_1);
  tree->Branch("event_2",&event_2);
  tree->Branch("overlap_12",&overlap_12);
  tree->Branch("overlap_12_norm",&overlap_12_norm);
  tree->Branch("wais_density_1",&wais_density_1);
  tree->Branch("wais_density_2",&wais_density_2);
  tree->Branch("wais_distance_1",&wais_distance_1);
  tree->Branch("wais_distance_2",&wais_distance_2);
  tree->Branch("wais_density_1_norm",&wais_density_1_norm);
  tree->Branch("wais_density_2_norm",&wais_density_2_norm);
  tree->Branch("wais_distance_1_norm",&wais_distance_1_norm);
  tree->Branch("wais_distance_2_norm",&wais_distance_2_norm);
  tree->Branch("L",&L); 
  tree->Branch("snr_1",&snr_1);
  tree->Branch("snr_2",&snr_2);
  tree->Branch("dphi1",&dphi1);
  tree->Branch("dphi2",&dphi2);
  tree->Branch("dtheta1",&dtheta1);
  tree->Branch("dtheta2",&dtheta2);
  tree->Branch("sum_1",&sum_1);
  tree->Branch("sum_2",&sum_2);

  while (tree->GetEntries() < N) 
  {

    UCorrelator::PointingResolution pointA; 
    UCorrelator::PointingResolution pointB; 

    int first_entry = gRandom->Integer(Nwais); 
    c.GetEntry(first_entry); 

    if (fabs(FFTtools::wrap(sum->peak[0][0].phi - sum->wais.phi,360,0)) > 5) continue; 
    if (fabs(sum->peak[0][0].theta - sum->wais.theta) > 4) continue; 
 //   printf("1st: %d\n", first_entry); 

    UCorrelator::ProbabilityMap map (&p);
    map.add(sum, pat, AnitaPol::kHorizontal, 0); 

    sum_1 = get_sum(g->NSegments(), map.getProbSums(false)); 

    m.computePointingResolution(sum, AnitaPol::kHorizontal, 0, &pointA); 

    wais_density_1 = map.getBaseSums()[293]; 
    wais_density_1_norm = map.getBaseSums(true)[293]; 
    wais_distance_1 = UCorrelator::ProbabilityMap::dens2dist(wais_density_1, UCorrelator::ProbabilityMap::get_two_pi_sqrt_det(pointA.getdPhi(), pointA.getdTheta(), pointA.getCorr())); 
    wais_distance_1_norm = UCorrelator::ProbabilityMap::dens2dist(wais_density_1_norm, UCorrelator::ProbabilityMap::get_two_pi_sqrt_det(pointA.getdPhi(), pointA.getdTheta(), pointA.getCorr())); 

    snr_1 = sum->trainingCoherent().snr; 
    run_1 = sum->run; 
    event_1 = sum->eventNumber; 

    double phiAA = sum->peak[0][0].phi; 
    double thetaAA = sum->peak[0][0].theta; 

    UsefulAdu5Pat patA(pat); 

    double latA = sum->peak[0][0].latitude; 
    double lonA = sum->peak[0][0].longitude; 
    double altA = sum->peak[0][0].altitude; 

    int second_entry = gRandom->Integer(Nwais); 
    if (first_entry == second_entry) continue; 

    c.GetEntry(second_entry); 

    UsefulAdu5Pat patB(pat); 

    double phiBB = sum->peak[0][0].phi; 
    double thetaBB = sum->peak[0][0].theta; 
    double latB = sum->peak[0][0].latitude; 
    double lonB = sum->peak[0][0].longitude; 
    double altB = sum->peak[0][0].altitude; 

    snr_2 = sum->trainingCoherent().snr; 
    run_2 = sum->run; 
    event_2 = sum->eventNumber; 
    if (fabs(FFTtools::wrap(sum->peak[0][0].phi - sum->wais.phi,360,0)) > 5) continue; 
    if (fabs(sum->peak[0][0].theta - sum->wais.theta) > 4) continue; 
//    printf("2nd: %d\n", second_entry); 

    map.add(sum, pat, AnitaPol::kHorizontal, 0); 
    sum_2 = get_sum(g->NSegments(), map.getProbSums(false)) - sum_1; 

    m.computePointingResolution(sum, AnitaPol::kHorizontal, sum->trainingPeakInd(), &pointB); 
    wais_density_2 = map.getBaseSums()[293] - wais_density_1; 
    wais_density_2_norm = map.getBaseSums(true)[293] - wais_density_1_norm; 
    wais_distance_2 = UCorrelator::ProbabilityMap::dens2dist(wais_density_2, UCorrelator::ProbabilityMap::get_two_pi_sqrt_det(pointB.getdPhi(), pointB.getdTheta(), pointB.getCorr())); 
    wais_distance_2_norm = UCorrelator::ProbabilityMap::dens2dist(wais_density_2_norm, UCorrelator::ProbabilityMap::get_two_pi_sqrt_det(pointB.getdPhi(), pointB.getdTheta(), pointB.getCorr())); 

    c.GetEntry(first_entry); 
    overlap_12 = map.overlap(sum, pat, AnitaPol::kHorizontal, 0,false);
    overlap_12_norm = map.overlap(sum, pat, AnitaPol::kHorizontal, 0, true); 


    double phiAB, phiBA, thetaAB, thetaBA ;
    patA.getThetaAndPhiWave(lonB, latB, altB, thetaAB, phiAB); 
    patB.getThetaAndPhiWave(lonA, latA, altA, thetaBA, phiBA); 
    phiAB *= TMath::RadToDeg();
    phiBA *= TMath::RadToDeg();
    thetaAB *= TMath::RadToDeg();
    thetaBA *= TMath::RadToDeg();

    dphi1 = pointA.getdPhi(); 
    dphi2 = pointB.getdPhi(); 
    dtheta1 = pointA.getdTheta(); 
    dtheta2 = pointB.getdTheta(); 

    L = sqrt(
        TMath::Power( FFTtools::wrap(phiAA-phiAB,360,0) / pointA.getdPhi(),2)  + 
        TMath::Power( FFTtools::wrap(phiBB-phiBA,360,0) / pointB.getdPhi(),2)  + 
        TMath::Power( (thetaAA-thetaAB) / pointA.getdTheta(),2)  + 
        TMath::Power( (thetaBB-thetaBA) / pointB.getdTheta(),2) );


    printf("evA: %d evB: %d (pair index %d)\n", event_1, event_2, tree->GetEntries()); 
    printf("  L=%g, overlap: %g, overlap_norm: %g\n", L, overlap_12, overlap_12_norm); 
    printf("    phiAA,phiAB, sigmaAphi: (%g %g, %g)  phiBA, phiBB, sigmaBphi:  (%g %g, %g)\n",phiAA, phiAB,pointA.getdPhi(), phiBA, phiBB, pointB.getdPhi()); 
    printf("    thetaAA,thetaAB, sigmaAtheta, (%g %g, %g)  thetaBA, thetaBB, sigmaBtheta(%g %g, %g)\n",thetaAA, thetaAB,pointA.getdTheta(), thetaBA, thetaBB, pointB.getdTheta()); 

    fout.cd(); 
    tree->Fill(); 
  }

  tree->Write(); 


  return 0; 
}
