
/** Example macro to grab waveforms.... It then writes them to a file. */ 

#include "AnitaVersion.h" 

void grabWaveforms(int run, int howmany = 10, const char * outfile = "example.root") 
{

  //May be needed if the default is not ANITA 3
  AnitaVersion::set(3); 

  //initialize the dataset 
  AnitaDataset d(run); 


  //initialize the analysis config
  UCorrelator::AnalysisConfig cfg; 

  //analyze hpol only... (default is to analyze both pols) 
  cfg.start_pol = AnitaPol::kHorizontal; 
  cfg.end_pol = AnitaPol::kHorizontal; 

  // Enable deconvolution with individual channel responses 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 

  // Use "allpass" deconvolution (not really a deconvolution) 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 


  //initialize the analyzer
  UCorrelator::Analyzer analyzer(&cfg,true); 


  //initialize the filter strategy
  FilterStrategy strategy; 
  UCorrelator::fillStrategyWithKey(&strategy,"sinsub_10_3_ad_2"); 



  //these will temporarily store the waveforms  and event number
  std::vector<int> eventNumbers; 
  std::vector<TGraph *> waveforms; 
  std::vector<TGraph *> waveforms_xpol; 

  for (int i = 0; i < d.N(); i++) 
  {

    //stop when we have enough
    if (waveforms.size() >= howmany) break; 
    
    //load the event from the dataset 
    d.getEntry(i); 

    //we will cut on WAIS pulses here, but you can use whateer condition you want
    UsefulAdu5Pat pat(d.gps()); 
    if (!UCorrelator::isWAISHPol(&pat, d.header())) continue; 
    

    printf("Using entry %d (event %d))", i, d.header()->eventNumber); 
    //create the filtered anita event
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

    //the summary, which we'll actually ignore here... 
    AnitaEventSummary sum; 

    //analyze it
    analyzer.analyze(&ev, &sum); 


    //save the deconvolved waveforms and event number
    eventNumbers.push_back(d.header()->eventNumber); 
    waveforms.push_back( new TGraph( *(analyzer.getDeconvolved(AnitaPol::kHorizontal,0)->even()) ) ); 
    waveforms_xpol.push_back( new TGraph (*(analyzer.getDeconvolvedXpol(AnitaPol::kHorizontal,0)->even()) ) ); 
  }



  //Now let's write these all to a file . This isn't the best way to do this, but this is just an example. 

  TFile out (outfile, "RECREATE"); 

  for (size_t i = 0; i < eventNumbers.size(); i++) 
  {
    //just to format the name
    TString name; 
    name.Form("g%d", eventNumbers[i]); 
    waveforms[i]->Write(name.Data()); 
    name.Form("g%d_xpol", eventNumbers[i]); 
    waveforms_xpol[i]->Write(name.Data()); 
     
  }

  

}

