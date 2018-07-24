#include "PolarityMachine.h"
#include "TimeDependentAverage.h"
#include "SystemResponse.h"
#include "AnitaTemplates.h"

  PolarityMachine::PolarityMachine()
: wfcomb(15, 3, true, false, 0)
{
  length = 2048;
  lengthFFT = length/2 + 1;
  dT = 0.1;
  dF = 1./(dT*length);
  kTemplatesLoaded = false;
  windowSizeNs = 10;

  d = new AnitaDataset(200);
  deconv = new AnitaResponse::AllPassDeconvolution;
  fNotchStr = "";

  for (int i=0; i < numCRTemplates; i++) {
    simulatedCRs[i] = NULL;
  }
  if(AnitaVersion::get() == 4) fillNotchConfigs();

  zeroInternals();
}

PolarityMachine::~PolarityMachine() {
  zeroInternals();
}


void PolarityMachine::zeroInternals() {
  /* in case you want to delete everything */

  if(kTemplatesLoaded)
  {
    for (int i=0; i < numCRTemplates; i++) {
      delete simulatedCRs[i];
      simulatedCRs[i] = NULL;
    }
    delete theImpTemplate;
    delete CRtemplate;
    delete CRtemplate_deconv;
    delete CRtemplate_deconv_windowed;
  } 

  kTemplatesLoaded = false;
  return;
}

void PolarityMachine::fillNotchConfigs()
{
  int tempTime;
  std::string tempStr;
  char* installDir = getenv("ANITA_UTIL_INSTALL_DIR");
  std::ifstream inf(Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/index.txt", installDir));
  while(inf >> tempStr >> tempTime)
  {
    payloadTimes.push_back(tempTime);
    notchConfigs.push_back(tempStr);
  }
  inf.close();
} 

void PolarityMachine::getCRTemplates(int version) {

  char* installDir = getenv("ANITA_UTIL_INSTALL_DIR");
  TString fname;
  TString CRtemplateName;
  TString CRtemplateName_deconv;
  if(version == 4)
  {
    fname = Form("%s/share/AnitaAnalysisFramework/templates/crTmpltsA4_%s.root", installDir, fNotchStr.c_str());
    CRtemplateName = Form("%s/share/AnitaAnalysisFramework/templates/averageA4CR_%s.txt", installDir, fNotchStr.c_str());
    CRtemplateName_deconv = Form("%s/share/AnitaAnalysisFramework/templates/deconvolvedCRA4average.txt", installDir);
  }

  else
  {
    fname = Form("%s/share/AnitaAnalysisFramework/templates/crTmpltsA3.root", installDir);
    CRtemplateName = Form("%s/share/AnitaAnalysisFramework/templates/averageA3CR.txt", installDir);
    CRtemplateName_deconv = Form("%s/share/AnitaAnalysisFramework/templates/deconvolvedCRA3average.txt", installDir);
  }

  TFile *inFile = TFile::Open(fname.Data());
  std::stringstream name;

  for (int i=0; i<numCRTemplates; i++) {
    //want to get graphs 13 through 24 (like in makeTemplate.C)
    int wave = i+13; //peak seems to be at around the 13th one, then by 23 it is basically zero
    name.str("");
    name << "disp" << wave;
    TGraph* guncut = (TGraph*)inFile->Get(name.str().c_str());
    TGraph* g = WindowingTools::windowCut(guncut, length);
    g = FFTtools::normalizeWaveform(g);
    delete guncut;
    if(version == 4) //the tuffs cause a flip of 180 vs the a3 stuff
    {
      for(int j = 0; j < g->GetN(); j++) g->GetY()[j] *= -1;
    }

    simulatedCRs[i] = new AnalysisWaveform(g->GetN(), g->GetX(), g->GetY(), dT);
    simulatedCRPeaks[i] = FFTtools::getPeakVal((TGraph*) simulatedCRs[i]->hilbertEnvelope());
    delete g;
  }
  TGraph* gCR = new TGraph(CRtemplateName.Data());
  gCR = FFTtools::normalizeWaveform(gCR);
  CRtemplate = new AnalysisWaveform(gCR->GetN(), gCR->GetX(), gCR->GetY(), dT);

  TGraph* gCR_deconv = new TGraph(CRtemplateName_deconv.Data());
  gCR_deconv = FFTtools::normalizeWaveform(gCR_deconv);
  CRtemplate_deconv = new AnalysisWaveform(gCR_deconv->GetN(), gCR_deconv->GetX(), gCR_deconv->GetY(), dT);
  int temp = -1;
  //gCR_deconv = WindowingTools::windowWave(gCR_deconv, temp, windowStart, windowStart, windowEnd, windowEnd);
  delete gCR_deconv;
  gCR_deconv = windowWaveform(CRtemplate_deconv, windowSizeNs);
  gCR_deconv = FFTtools::normalizeWaveform(gCR_deconv);
  CRtemplate_deconv_windowed = new AnalysisWaveform(gCR_deconv->GetN(), gCR_deconv->GetX(), gCR_deconv->GetY(), dT);

  delete gCR;
  delete gCR_deconv;

  inFile->Close();

  return;
}

void PolarityMachine::getImpulseResponseTemplate(int version) {

  char* installDir = getenv("ANITA_UTIL_INSTALL_DIR");
  TString name;
  if(version == 4)
  {
    name = Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/averages/%s.imp", installDir, fNotchStr.c_str());
  }
  else name = Form("%s/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp", installDir);
  TGraph *gunpad = new TGraph(name.Data());
  TGraph *g = FFTtools::padWaveToLength(gunpad, length);

  theImpTemplate = new AnalysisWaveform(g->GetN(), g->GetX(), g->GetY(), dT);
  delete g;
  delete gunpad;
  return;
}


void PolarityMachine::loadTemplates(unsigned int eventTime, int version) {

  std::cout << "Loading templates, length=" << length << std::endl;
  std::string tempStr = "";
  int i = 0;
  if(version == 4)
  {
    while(i < payloadTimes.size())
    {
      if(eventTime < payloadTimes.at(i)) break;
      i++;
    }
    tempStr = notchConfigs.at(i);
    if(kTemplatesLoaded)
    {
      if(tempStr.compare(fNotchStr) == 0) return;
    }
  }
  fNotchStr = tempStr;
  if (kTemplatesLoaded) zeroInternals();
  getCRTemplates(version);
  getImpulseResponseTemplate(version);

  kTemplatesLoaded = true;
  return;
}

TH2D* PolarityMachine::generatePolarityMeasurements(int N, int eventNumber, double cr_phi, double cr_theta, double cr_snr, int pol, int metric,  const char * outhistoname,  std::vector<int> whichTemplates)
{
  d->getEvent(eventNumber, true);
  //make sure templates are loaded
  loadTemplates(d->header()->realTime);
  if(cr_snr <= 0)
  {
    printf("something is wrong with your SNR !!\n");
    return 0;
  }

  TRandom3* tr = new TRandom3(0);
  int j = 0;
  TH2D* h = new TH2D(outhistoname, outhistoname, 220, -2.1,2.1, 220, -2.1, 2.1);
  h->GetXaxis()->SetTitle("coherentPolarity");
  h->GetYaxis()->SetTitle("deconvolvedPolarity");

  double average_coh_pol = 0;
  double average_deco_pol = 0;
  for(int i = 0; i < N; i++)
  {
    //generate a noisy wf and test its polarity
    AnalysisWaveform* noiseWf = makeNoiseWaveformFromMinBias(eventNumber, cr_phi, cr_theta, pol, i);
    double rms = TMath::RMS(noiseWf->even()->GetN(), noiseWf->even()->GetY());
    double wf_snr = simulatedCRPeaks[whichTemplates[j]] / rms;

    AnalysisWaveform* noisyWf = generateNoisyWaveform(simulatedCRs[whichTemplates[j]], noiseWf, cr_snr/wf_snr, tr);
    double coherent_polarity = testPolarity(metric, noisyWf, 0);

    //now deconvolve it and test the polarity of that
    deconv->deconvolve(noisyWf->Nfreq(), noisyWf->deltaF(), noisyWf->updateFreq(), theImpTemplate->freq());
    double deco_polarity = testPolarity(metric, noisyWf, 1);

    average_coh_pol += coherent_polarity;
    average_deco_pol += deco_polarity;

    h->Fill(coherent_polarity, deco_polarity);
    j = (j < whichTemplates.size()-1) ? j+1 : 0; //cycle through the templates you want to look at
    delete noiseWf;
    delete noisyWf;
  }
  printf("average coherent polarity = %g, average deconvolved polarity  = %g\n", average_coh_pol/N, average_deco_pol/N);

  return h;
}

AnalysisWaveform* PolarityMachine::generateNoisyWaveform(AnalysisWaveform* templateWf, AnalysisWaveform* noiseWf, double snr, TRandom3* tr)
{
  double new_snr = tr->Gaus(snr, 1./snr);
  AnalysisWaveform* outWf = new AnalysisWaveform(templateWf->even()->GetN(), templateWf->even()->GetY(), templateWf->deltaT(), 0);

  for(int i = 0; i < noiseWf->Neven(); i++)
    outWf->updateEven()->GetY()[i] = templateWf->even()->GetY()[i] * snr + noiseWf->even()->GetY()[i];

  return outWf;
}

double PolarityMachine::testPolarity(int metric, AnalysisWaveform* wf, bool deconvolved)
{
  TGraph* gwf = gwf = new TGraph(wf->even()->GetN(), wf->even()->GetX(), wf->even()->GetY());
  TGraph* gCRwf = 0;
  TGraph* gCorr = 0;
  if(metric == 4)
  {
    double maxVal = TMath::MaxElement(gwf->GetN(), gwf->GetY());
    double minVal = TMath::MinElement(gwf->GetN(), gwf->GetY());
    double retVal = (maxVal + minVal) / (maxVal - minVal);
    delete gwf;
    return retVal;
  }
  else if(!deconvolved || metric < 2)
  {
    gCRwf = deconvolved ? new TGraph(CRtemplate_deconv->even()->GetN(), CRtemplate_deconv->even()->GetX(), CRtemplate_deconv->even()->GetY()) : new TGraph(CRtemplate->even()->GetN(), CRtemplate->even()->GetX(), CRtemplate->even()->GetY());
  }
  else if(metric == 2 || metric == 3)
  {
    gCRwf = new TGraph(CRtemplate_deconv_windowed->even()->GetN(), CRtemplate_deconv_windowed->even()->GetX(), CRtemplate_deconv_windowed->even()->GetY());
    int temp = -1;
    //windowing could probably be tuned up
    //gwf = WindowingTools::windowWave(gwf, temp, windowStart, windowStart, windowEnd, windowEnd);
    delete gwf;
    gwf = windowWaveform(wf, windowSizeNs);
  }
  gCorr = FFTtools::getCorrelationGraph(gwf, gCRwf);
  double normalization = 1./(gwf->GetRMS(2) * gCRwf->GetRMS(2) * gwf->GetN()/TMath::Power(2, int(TMath::Log2(gwf->GetN()))));
  double maxCorr = normalization * TMath::MaxElement(gCorr->GetN(), gCorr->GetY());
  double minCorr = normalization * TMath::MinElement(gCorr->GetN(), gCorr->GetY());
  double peak2sidelobe = abs(maxCorr/minCorr);

  delete gwf;
  delete gCRwf;
  delete gCorr;

  if(metric%2 == 1) return peak2sidelobe;
  if(abs(maxCorr) > abs(minCorr)) return 1;
  return -1;
}

AnalysisWaveform* PolarityMachine::makeNoiseWaveformFromMinBias(int eventNumber, double phi, double theta, int pol, int current_N)
{
  //TODO make the output of this be 2048 points long or 2046 in awf form or change everything else to 1500 points
  d->getEvent(eventNumber, true);
  if(current_N%2 == 0)
  {
    for(int i = 0; i < (current_N/2) + 1; i++)
      d->nextMinBiasEvent();
  }
  if(current_N%2 == 1)
  {
    for(int i = 0; i < (current_N/2) + 1; i++)
      d->previousMinBiasEvent();
  }
  FilterStrategy strat;

  FilteredAnitaEvent fae(d->useful(), &strat, d->gps(), d->header());
  wfcomb.combine(phi, theta, &fae, AnitaPol::AnitaPol_t(pol), 0, -25, 125);
  TGraph* g = new TGraph(wfcomb.getCoherent()->even()->GetN(), wfcomb.getCoherent()->even()->GetX(), wfcomb.getCoherent()->even()->GetY());
  g = FFTtools::normalizeWaveform(g);
  AnalysisWaveform* theOut = new AnalysisWaveform(g->GetN(), g->GetX(), g->GetY(), dT);
  delete g;

  return theOut;
}

TGraph* PolarityMachine::windowWaveform(AnalysisWaveform* wf, double window_size_ns)
{
  double * x = wf->even()->GetY();
  double * xh = wf->hilbertTransform()->even()->GetY();
  double max_dI = 0;
  int max_ind = -1;
  for(int i = 0; i < wf->Neven(); i++)
  {
    std::complex<double> X(x[i], xh[i]);
    double dI = std::norm(X);
    if(dI > max_dI)
    {
      max_dI = dI;
      max_ind = i;
    }
  }
  int window_edge = int(ceil(window_size_ns/(2*wf->deltaT())));
  TGraph* theOut = new TGraph(wf->even()->GetN(), wf->even()->GetX(), wf->even()->GetY());
  double alpha = 2 * (1. - double(2 * window_edge)/double(theOut->GetN())); //tukey window alpha
  for(int i = 0; i < theOut->GetN(); i++)
  {
    if(i >= (max_ind - window_edge) && i <= (max_ind + window_edge)) continue;
    else if (i < (max_ind - window_edge))
    {
      double tukey_const = double(2*i)/(alpha * double(theOut->GetN()-1));
      double coeff = .5*(1 + TMath::Cos(TMath::Pi()*(tukey_const-1)));
      theOut->GetY()[i] *= coeff;
    }
    else
    {
      int j = 2 * max_ind - i;
      double tukey_const = double(2*j)/(alpha * double(theOut->GetN()-1));
      double coeff = .5*(1 + TMath::Cos(TMath::Pi()*(tukey_const-1)));
      theOut->GetY()[i] *= coeff;
    }
  }
  return theOut;
}




