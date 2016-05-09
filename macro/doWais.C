
void doWais(int run = 352, int max = 0, bool sine_subtract = false)
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  cfg.start_pol = AnitaPol::kHorizontal; 
  cfg.end_pol = AnitaPol::kHorizontal; 
  

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 
  if (max) outname.Form("wais/wais_hpol_%d_max_%d%s.root",run,max, sine_subtract ? "_sinsub" : "" ); 
  else outname.Form("wais/wais_hpol_%d%s.root",run, sine_subtract ? "_sinsub" : "" ); 

  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree("wais","WAIS Hpol"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 


  FilterStrategy strategy (&ofile); 
  if (sine_subtract) 
  {
//    double fmins[2] = {0.23, 0.43}; 
//    double fmaxs[2] = {0.29, 0.49}; 
    strategy.addOperation(new UCorrelator::SineSubtractFilter(0.05, 0, 4)); 
  }
  else
  {
    UCorrelator::applyAbbysFilterStrategy(&strategy); 
  }

//  printf("Strategy applied!\n"); 

  double waisPhiExpected, waisThetaExpected; 
  tree->Branch("summary",&sum); 
  tree->Branch("waisPhiExpected", &waisPhiExpected); 
  tree->Branch("waisThetaExpected", &waisThetaExpected); 

  int ndone = 0; 
  for (int i =0 ; i < d.N(); i++)
  {

    d.getEntry(i); 
    printf("----(%d)-----\n",i); 
    UsefulAdu5Pat pat(d.gps()); 

    if (UCorrelator::isWAISHPol(&pat, d.header()))
    {
      pat.getThetaAndPhiWaveWaisDivide(waisThetaExpected, waisPhiExpected); 
      waisThetaExpected *= 180 / TMath::Pi(); 
      waisPhiExpected *= 180 / TMath::Pi(); 
      printf("Processing event %d (%d)\n",d.header()->eventNumber,ndone); 
      FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

      analyzer.analyze(&ev, sum); 
      ofile.cd(); 
      tree->Fill(); 
      ndone++; 
    }

    if (max && ndone > max) break; 

  }

  ofile.cd(); 
  tree->Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
}
