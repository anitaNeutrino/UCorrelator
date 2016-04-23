
void doWais(int run = 352)
{

  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  cfg.start_pol = AnitaPol::kHorizontal; 
  cfg.end_pol = AnitaPol::kHorizontal; 

  UCorrelator::Analyzer analyzer(&cfg); 

  TFile ofile(TString::Format("wais/wais_hpol_%d.root", run), "RECREATE"); 
  TTree * tree = new TTree("wais","WAIS Hpol"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 


  FilterStrategy strategy (&ofile); 
  UCorrelator::applyAbbysFilterStrategy(&strategy); 

  printf("Strategy applied!\n"); 

  double waisPhiExpected, waisThetaExpected; 
  tree->Branch("summary",&sum); 
  tree->Branch("waisPhiExpected", &waisPhiExpected); 
  tree->Branch("waisThetaExpected", &waisThetaExpected); 

  int ndone = 0; 
  for (int i =0 ; i < d.N(); i++)
  {

    d.getEntry(i); 
//    printf("%d\n",i); 
    UsefulAdu5Pat pat(d.gps()); 

    if (UCorrelator::isWAISHPol(&pat, d.header()))
    {
      pat.getThetaAndPhiWaveWaisDivide(waisThetaExpected, waisPhiExpected); 
      printf("Processing event %d (%d)\n",d.header()->eventNumber,ndone); 
      FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

      analyzer.analyze(&ev, sum); 
      ofile.cd(); 
      tree->Fill(); 
      ndone++; 
    }

    if (ndone > 50) break; 

  }

  ofile.cd(); 
  tree->Write(); 

}
