
void doDecimated(int run = 352, int max = 0)
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run,true); 
  UCorrelator::AnalysisConfig cfg; 
  

  UCorrelator::Analyzer analyzer(&cfg); 

  TFile ofile(TString::Format("decimated/%d.root", run), "RECREATE"); 
  TTree * tree = new TTree("decimated","decimated"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 


  FilterStrategy strategy (&ofile); 
  UCorrelator::applyAbbysFilterStrategy(&strategy); 

//  printf("Strategy applied!\n"); 

  double waisPhiExpected, waisThetaExpected; 
  tree->Branch("summary",&sum); 
  tree->Branch("waisPhiExpected", &waisPhiExpected); 
  tree->Branch("waisThetaExpected", &waisThetaExpected); 

  int ndone = 0; 
  for (int i =0 ; i < d.N(); i++)
  {

    d.getEntry(i); 
//    printf("%d\n",i); 
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 
    analyzer.analyze(&ev, sum); 
    ofile.cd(); 
    tree->Fill(); 
    ndone++; 

    if (max && ndone > max) break; 

  }

  ofile.cd(); 
  tree->Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
}
