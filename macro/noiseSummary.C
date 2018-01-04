void noiseSummary(int run=300, int numEntries=1000) {

  cout << " Hello! :D" << endl;

  //open data
  AnitaDataset *data = new AnitaDataset(run);

  //determine number of events
  int lenEntries = data->N();
  if (numEntries == -1 || numEntries > lenEntries) {
    numEntries = lenEntries;
  }
  cout << "Found " << lenEntries << " events in run " << run << ".  Using " << numEntries << " of them." << endl;

  stringstream name;


  //open output files
  //for the event summary
  string outFileName = "eventSummary";
  name.str(""); 
  name << outFileName << run << ".root";
  cout << "Making output file " << name.str() << " and output tree and filling it with stuff" << endl;
  TFile *outFile = TFile::Open(name.str().c_str(),"recreate");
  //make output tree
  TTree* summaryTree = new TTree("summaryTree","summaryTree");
  //make output things
  AnitaEventSummary *eventSummary = new AnitaEventSummary();
  summaryTree->Branch("eventSummary",&eventSummary);

  //and the noise summary (you can combine them if you want)
  outFileName = "noiseSummary";
  name.str(""); 
  name << outFileName << run << ".root";
  cout << "Making output file " << name.str() << " and output tree and filling it with stuff" << endl;
  TFile* noiseOutFile = TFile::Open(name.str().c_str(),"recreate");
  TTree *noiseTree = new TTree("noiseTree","noiseTree");
  AnitaNoiseSummary *noiseSummary = new AnitaNoiseSummary();
  noiseTree->Branch("noiseSummary",&noiseSummary);

  //make the noise machine
  AnitaNoiseMachine *noiseMachine = new AnitaNoiseMachine();
  noiseMachine->fillMap = true;
  
  //make the analyzer  
  cout << "Making the Analyzer" << endl;
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 
  config->response_option = UCorrelator::AnalysisConfig::ResponseOption_t::ResponseSingleBRotter;
  //and create an analyzer object.  Currently the "true" is required for the noiseMachine (interactiveMode)
  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true);

  
  //make a filtering plan
  FilterStrategy *strategy = new FilterStrategy();
  char* specAvgDir = getenv("UCORRELATOR_SPECAVG_DIR");
  const UCorrelator::TimeDependentAverageLoader *specAvgLoader = new UCorrelator::TimeDependentAverageLoader(specAvgDir);
  UCorrelator::AdaptiveBrickWallFilter *brickWall = new UCorrelator::AdaptiveBrickWallFilter(specAvgLoader,2,false);
  strategy->addOperation(brickWall);


  //loop through events!
  cout << "Okay lets go" << endl;
  for (int entry=0; entry<numEntries; entry++) {
    data->getEntry(entry);

    //only want minbias stuff, skip it otherwise
    int trigType = data->header()->trigType&0x0F;
    if (trigType == 1) continue;

    //say hello and how far you are
    cout << entry << " / " << numEntries << endl;
      
    //make your filtered event
    FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data->useful(), strategy, data->gps(), data->header());

    //clear the eventSummary so that I can fill it up with the analyzer
    eventSummary->zeroInternals();

    //then analyze the filtered event (analyzer fills up most of eventSummary too)!
    analyzer->analyze(filteredEvent, eventSummary); 

    //do the noise analysis (only update it if it is a min bias though 
    //unless you want to see what biased triggered noise looks like
    if (trigType == 1) {
      noiseSummary->isMinBias = false;
    }
    else {
      noiseSummary->isMinBias = true;
      noiseMachine->updateMachine(analyzer,filteredEvent);
    }

    //use the noiseMachine to fill up the noiseSummary and last few pieces of eventSummary
    noiseMachine->fillNoiseSummary(noiseSummary);
    noiseMachine->fillEventSummary(eventSummary);

    //fill the trees!
    summaryTree->Fill();
    noiseTree->Fill();

    //delete the filteredEvent from memory or you'll end up leaking it everywhere
    delete filteredEvent;
    //You're in interactive mode, so you have to clear the memory from the analyzer
    analyzer->clearInteractiveMemory();

    
  }

  //write everything out and close the files
  outFile->cd();
  summaryTree->Write();
  outFile->Close();

  noiseOutFile->cd();
  noiseTree->Write();
  noiseOutFile->Close();


  //bye!
  cout << "Done!  Thanks :) " << endl;


  return ;
}
