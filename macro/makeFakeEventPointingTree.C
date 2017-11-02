#include "FFTtools.h" 

void makeFakeEventPointingTree(int run = 342, bool replace_noise = false, bool make_images = false) 
{

 FFTtools::loadWisdom("wisdom.dat"); 
 
  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 3; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
  cfg.use_forced_trigger_rms = !replace_noise; 


  AnitaDataset d(run); 

  UCorrelator::Analyzer  analyzer(&cfg,make_images); 

  system("mkdir -p fake"); 

  if (make_images) 
  {
    TString str; 
    str.Form("mkdir -p fake/%d/", run); 
    system(str.Data()); 
  }


  TFile fout(TString::Format("fake/fakeEvents%d_%d.root",run, replace_noise),"RECREATE"); 

  AnitaEventFaker faker("IndividualBRotter"); 
  TTree * tree = new TTree("fake","Fake events"); 

  double truePhi, trueTheta; 
  double JonesH, JonesV; 
  double A; 
  AnitaEventSummary * sum = new AnitaEventSummary; 


  FilterStrategy strategy; 
  if (!replace_noise) UCorrelator::fillStrategyWithKey(&strategy, "sinsub_10_3_ad_2"); 

  TruthAnitaEvent truth; 

  tree->Branch("summary",&sum); 
  tree->Branch("truePhi",&truePhi); 
  tree->Branch("trueTheta",&trueTheta); 
  tree->Branch("JonesH",&JonesH); 
  tree->Branch("JonesV",&JonesV); 
  tree->Branch("A",&A); 

  TCanvas * ca = 0;

  if (make_images) 
  {
    ca = new TCanvas("ca","ca",1920,1200); 
    ca->Divide(1,2); 
  }

  for (int i = 0; i < d.N(); i++) 
  {
    d.getEntry(i); 
    if (d.header()->trigType & 1) continue; //only force triggers  
    printf("%d\n",i); 
    if (replace_noise) faker.makePureNoiseEvent(0.1, d.useful()); 
    A = gRandom->Uniform(1,30); 
    truePhi = gRandom->Uniform(0,360); 
    trueTheta = gRandom->Uniform(-20,50); 
    truth.payloadPhi = truePhi; 
    truth.payloadTheta = trueTheta; 

    JonesH = gRandom->Uniform(0,1); 
    JonesV = gRandom->Uniform(0,1); 

    double mag = sqrt(JonesH*JonesH + JonesV*JonesV); 
    JonesH/= mag; 
    JonesV/= mag; 
    faker.addSignal(d.useful(), trueTheta, truePhi, A, std::complex<double>(JonesH,0), std::complex<double>(JonesV,0)); 


    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 
    analyzer.analyze(&ev, sum, &truth);
    if (make_images)
    {
      analyzer.drawSummary((TPad*) ca->GetPad(1), (TPad*) ca->GetPad(2)); 
      ca->SaveAs(TString::Format("fake/%d/%d.png",run, replace_noise ? (int) tree->GetEntries() : d.header()->eventNumber) );
    }
    fout.cd(); 
    tree->Fill(); 

  }

  fout.cd(); 
  tree->Write(); 
}

