#include "PolarityMachine.h"
#include "TimeDependentAverage.h"
#include "SystemResponse.h"
#include "AnitaTemplates.h"

PolarityMachine::PolarityMachine(int padX)
: wfcomb(15, 3, true, false, 0)
{
  length = 2048;
  lengthFFT = length/2 + 1;
  dT = 0.1;
  dF = 1./(dT*length);
  kTemplatesLoaded = false;
  windowSizeNs = 10;
  padFactor = padX;
  test_coherent = true;
  forwardMinBias = 0;
  previousMinBias = 0;

  d = new AnitaDataset(41);
  deconv = new AnitaResponse::AllPassDeconvolution;
  fNotchStr = "";
  inputWaveform = new TGraph;
  crWaveform = new TGraph;
  correlationGraph = new TGraph;

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
  forwardMinBias = 0;
  previousMinBias = 0;

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
  CRtemplate->padFreq(padFactor);

  TGraph* gCR_deconv = new TGraph(CRtemplateName_deconv.Data());
  gCR_deconv = FFTtools::normalizeWaveform(gCR_deconv);
  CRtemplate_deconv = new AnalysisWaveform(gCR_deconv->GetN(), gCR_deconv->GetX(), gCR_deconv->GetY(), dT);
  CRtemplate_deconv->padFreq(padFactor);
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

  for(int i = 0; i < N; i++)
  {
    //generate a noisy wf and test its polarity
    AnalysisWaveform* noiseWf = makeNoiseWaveformFromMinBias(eventNumber, cr_phi, cr_theta, pol, i);
    double rms = TMath::RMS(noiseWf->even()->GetN(), noiseWf->even()->GetY());
    double wf_snr = simulatedCRPeaks[whichTemplates[j]] / rms;

    AnalysisWaveform* noisyWf = generateNoisyWaveform(simulatedCRs[whichTemplates[j]], noiseWf, cr_snr/wf_snr, tr);
    double coherent_polarity = (test_coherent) ? testPolarity(metric, noisyWf, 0) : 0;

    //now deconvolve it and test the polarity of that
    deconv->deconvolve(noisyWf->Nfreq(), noisyWf->deltaF(), noisyWf->updateFreq(), theImpTemplate->freq());
    double deco_polarity = testPolarity(metric, noisyWf, 1);

    h->Fill(coherent_polarity, deco_polarity);
    j = (j < whichTemplates.size()-1) ? j+1 : 0; //cycle through the templates you want to look at
    delete noiseWf;
    delete noisyWf;
  }

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
  if(wf->deltaT() != dT) 
  {
    AnalysisWaveform wfTemp(wf->Neven(), wf->even()->GetX(), wf->even()->GetY(), dT);
    //delete wf;
    wf = new AnalysisWaveform(wfTemp.Neven(), wfTemp.even()->GetX(), wfTemp.even()->GetY(), dT);
  }
  wf->padFreq(padFactor);
  if(deconvolved)
  {
    if(CRtemplate_deconv->Neven() != wf->Neven()) makeSameSize(CRtemplate_deconv, wf);
  }
  else
  {
    if(CRtemplate->Neven() != wf->Neven()) makeSameSize(CRtemplate, wf);
  }
  TGraph* gwf = new TGraph(wf->even()->GetN(), wf->even()->GetX(), wf->even()->GetY());
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
  else if(deconvolved && metric == 5)
  {
    gCRwf = peterWindow(CRtemplate_deconv, windowSizeNs);
    delete gwf;
    gwf = peterWindow(wf, windowSizeNs);
    double maxCRnum = TMath::MaxElement(gCRwf->GetN(), gCRwf->GetY());
    double maxWFnum = TMath::MaxElement(gwf->GetN(), gwf->GetY());
  }
  else if(deconvolved && metric == 6)
  {
    delete gwf;
    gwf = peterWindow(wf, windowSizeNs, true);
    delete inputWaveform;
    inputWaveform = new TGraph(gwf->GetN(), gwf->GetX(), gwf->GetY());
    double ave_phase = 0;
    double n_samp = 0;
    for(int i = 0; i < gwf->GetN(); i++)
    {
      if(gwf->GetX()[i] > 0.6) break;
      if(gwf->GetX()[i] < 0.18) continue;
      n_samp++;
      ave_phase += gwf->GetY()[i];
    }
    ave_phase/=n_samp;
    ave_phase = (ave_phase > 0) ? 1 : -1;
    delete gwf;
    return ave_phase;
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
  double normalization = TMath::Power(2, int(TMath::Log2(gwf->GetN())))/(gwf->GetRMS(2) * gCRwf->GetRMS(2) * gwf->GetN());
  for(int i = 0; i < gCorr->GetN(); i++) gCorr->GetY()[i] *= normalization;
  double maxCorr = TMath::MaxElement(gCorr->GetN(), gCorr->GetY());
  double minCorr = TMath::MinElement(gCorr->GetN(), gCorr->GetY());
  double peak2sidelobe = abs(maxCorr/minCorr);
  double fixedMetric = gCorr->GetY()[gCorr->GetN()/2];
  //printf("%g, %g, %g, %g, %g\n", normalization, maxCorr, minCorr, peak2sidelobe, fixedMetric);

  delete inputWaveform;
  delete crWaveform;
  delete correlationGraph;
  inputWaveform = new TGraph(gwf->GetN(), gwf->GetX(), gwf->GetY());
  crWaveform = new TGraph(gCRwf->GetN(), gCRwf->GetX(), gCRwf->GetY());
  correlationGraph = new TGraph(gCorr->GetN(), gCorr->GetX(), gCorr->GetY());
  delete gwf;
  delete gCRwf;
  delete gCorr;

  if(deconvolved && metric == 5) return fixedMetric;
  if(metric%2 == 1) return peak2sidelobe;
  if(abs(maxCorr) > abs(minCorr)) return 1;
  return -1;
}

double PolarityMachine::fourierPhasePolarity(int window_type, AnalysisWaveform* wf)
{
  TGraph* theOut;
  if(window_type == 0) theOut = peterWindow(wf, windowSizeNs, true);
  else if(window_type == 1) theOut = hilbertWindow(wf, windowSizeNs, true, true, 0);
  else if(window_type == 2) theOut = hilbertWindow(wf, windowSizeNs, true, true, 1);
  else if(window_type == 3) theOut = hilbertWindow(wf, windowSizeNs, true, false);

  delete inputWaveform;
  inputWaveform = new TGraph(theOut->GetN(), theOut->GetX(), theOut->GetY());
  double ret = 0;
  double n_samp = 0;
  for(int i = 0; i < theOut->GetN(); i++)
  {
    if(theOut->GetX()[i] > 0.6) break;
    if(theOut->GetX()[i] < 0.18) continue;
    ret += theOut->GetY()[i];
    n_samp++;
  }
  ret /= n_samp;
  delete theOut;
  return ret;
}

AnalysisWaveform* PolarityMachine::makeNoiseWaveformFromMinBias(int eventNumber, double phi, double theta, int pol, int current_N)
{
  //TODO make the output of this be 2048 points long or 2046 in awf form or change everything else to 1500 points
  d->getEvent(eventNumber, true);
  if(current_N%2 == 0) 
  {
    if(forwardMinBias != 0)
    {
      d->getEvent(forwardMinBias);
    }
    else 
    {
      for(int i = 0; i < offset/2; i++) d->nextMinBiasEvent();
    }
    d->nextMinBiasEvent();
    forwardMinBias = d->header()->eventNumber;
  }
  else
  {
    if(previousMinBias != 0)
    {
      d->getEvent(previousMinBias);
    }
    else 
    {
      for(int i = 0; i < offset/2; i++) d->previousMinBiasEvent();
    }
    d->previousMinBiasEvent();
    previousMinBias = d->header()->eventNumber;
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
  double prev = x[max_ind];
  int zc_ind = -1;
  int zc = 0;
  for(int i = max_ind; i<wf->Neven(); i++)
  {
    std::complex<double> X(x[i], xh[i]);
    double dI = std::norm(X);
    if(prev < 0 && x[i] > 0) zc = 1;
    if(prev > 0 && x[i] < 0) zc = 1;
    if(zc)
    {
      max_dI = dI;
      zc_ind = i;
      break;
    }
    prev = x[i];
  }
  prev = x[max_ind];
  int sample_diff = zc_ind - max_ind;
  zc = 0;
  for(int i = max_ind; i > 0; i--)
  {
    std::complex<double> X(x[i], xh[i]);
    double dI = std::norm(X);
    if(prev < 0 && x[i] > 0) zc = 1;
    if(prev > 0 && x[i] < 0) zc = 1;
    if(zc)
    {
      if(max_ind - i < sample_diff) zc_ind = i;
      break;
    }
    prev = x[i];
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

TGraph* PolarityMachine::peterWindow(AnalysisWaveform* wf, double window_size_ns, bool roll)
{
  if(window_size_ns <= 0)
  {
    TGraph* theOut = new TGraph(wf->even()->GetN(), wf->even()->GetX(), wf->even()->GetY());
    return theOut;
  }

  double * x = wf->even()->GetY();
  double * xh = wf->hilbertTransform()->even()->GetY();
  double max_dI = 0;
  int max_ind = -1;
  double prev;

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
  prev = x[max_ind];
  int zc_ind = -1;
  int zc = 0;
  for(int i = max_ind; i<wf->Neven(); i++)
  {
    std::complex<double> X(x[i], xh[i]);
    double dI = std::norm(X);
    if(prev < 0 && x[i] > 0) zc = 1;
    if(prev > 0 && x[i] < 0) zc = 1;
    if(zc)
    {
      max_dI = dI;
      zc_ind = i;
      break;
    }
    prev = x[i];
  }
  prev = x[max_ind];
  int sample_diff = zc_ind - max_ind;
  zc = 0;
  for(int i = max_ind; i > 0; i--)
  {
    std::complex<double> X(x[i], xh[i]);
    double dI = std::norm(X);
    if(prev < 0 && x[i] > 0) zc = 1;
    if(prev > 0 && x[i] < 0) zc = 1;
    if(zc)
    {
      if(max_ind - i < sample_diff) zc_ind = i;
      break;
    }
    prev = x[i];
  }
  int window_edge = int(ceil(window_size_ns/(2*wf->deltaT())));
  TGraph* theOut = new TGraph(window_edge*2);
  for(int i = 0; i < wf->even()->GetN(); i++)
  {
    if(i < (zc_ind - window_edge) || i > (zc_ind + window_edge)) continue;
    else
    {
      if(roll) 
      {
        theOut->SetPoint(i-(zc_ind - window_edge), wf->even()->GetX()[i], wf->even()->GetY()[i]);
      }
      else
      {
        theOut->SetPoint(i-(zc_ind - window_edge), i-(zc_ind - window_edge), wf->even()->GetY()[i]);
      }
    }
  }

  if(roll)
  {
    TGraph* tempOut = new TGraph(theOut->GetN());
    for(int i = 0; i < theOut->GetN(); i++)
    {
      if(i < theOut->GetN()/2)
      {
        tempOut->SetPoint(i, theOut->GetX()[i], theOut->GetY()[i + theOut->GetN()/2]);
      }
      else 
      {
        tempOut->SetPoint(i, theOut->GetX()[i], theOut->GetY()[i - theOut->GetN()/2]);
      }
    }
    delete theOut;
    AnalysisWaveform* tawf = new AnalysisWaveform(tempOut->GetN(), tempOut->GetY(), tempOut->GetX()[1] - tempOut->GetX()[0], 0);
    //tawf->padEven(1);
    theOut = new TGraph(tawf->phase()->GetN(), tawf->phase()->GetX(), tawf->phase()->GetY());
    for(int i = 0; i < theOut->GetN(); i++)
    {
      theOut->GetY()[i] *= TMath::RadToDeg();
    }
  }

  return theOut;
}

TGraph* PolarityMachine::hilbertWindow(AnalysisWaveform* wf, double window_size_ns, bool roll, bool use_cdf, int how_to_use_cdf)
{
  if(window_size_ns <= 0)
  {
    TGraph* theOut = new TGraph(wf->even()->GetN(), wf->even()->GetX(), wf->even()->GetY());
    return theOut;
  }

  double * x = wf->even()->GetY();
  double * xh = wf->hilbertTransform()->even()->GetY();
  double max_dI = 0;
  int max_ind = -1;
  double prev;

  if(!use_cdf)
  {
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
  }
  else
  {
    TGraph* tempCDF = new TGraph(wf->Neven());
    double cdfmean = 0;
    for(int i = 0; i < wf->Neven(); i++)
    {
      std::complex<double> X(x[i], xh[i]);
      double dI = std::norm(X);
      if(i > 0) 
      {
        tempCDF->SetPoint(i, i, dI + tempCDF->GetY()[i-1]);
        cdfmean += (dI + tempCDF->GetY()[i-1])/double(wf->Neven());
      }
      else 
      {
        tempCDF->SetPoint(i, i, dI);
        cdfmean += dI/double(wf->Neven());
      }
    }
    
    if(how_to_use_cdf == 0)
    {
      for(int i = 0; i < tempCDF->GetN(); i++)
      {
        if(tempCDF->GetY()[i] > cdfmean)
        {
          max_ind = i;
          break;
        }
      }
      max_ind = (abs(tempCDF->GetY()[max_ind] - cdfmean) > abs(tempCDF->GetY()[max_ind-1] - cdfmean)) ? max_ind -1 : max_ind;
    }
    else if(how_to_use_cdf == 1)
    {
      double target = 0.9 * tempCDF->GetY()[tempCDF->GetN()-1];
      for(int i = tempCDF->GetN()-1; i > -1; i--)
      {
        if(tempCDF->GetY()[i] < target)
        {
          max_ind = i;
          break;
        }
      }
      max_ind = (abs(tempCDF->GetY()[max_ind] - target) > abs(tempCDF->GetY()[max_ind+1] - target)) ? max_ind -1 : max_ind;
    }
  }

  int window_edge = int(ceil(window_size_ns/(2*wf->deltaT())));
  TGraph* theOut = new TGraph(window_edge*2);
  for(int i = 0; i < wf->even()->GetN(); i++)
  {
    if(i < (max_ind - window_edge) || i > (max_ind + window_edge)) continue;
    else
    {
      if(roll) 
      {
        theOut->SetPoint(i-(max_ind - window_edge), wf->even()->GetX()[i], wf->even()->GetY()[i]);
      }
      else
      {
        theOut->SetPoint(i-(max_ind - window_edge), i-(max_ind-window_edge), wf->even()->GetY()[i]);
      }
    }
  }

  if(roll)
  {
    TGraph* tempOut = new TGraph(theOut->GetN());
    for(int i = 0; i < theOut->GetN(); i++)
    {
      if(i < theOut->GetN()/2)
      {
        tempOut->SetPoint(i, theOut->GetX()[i], theOut->GetY()[i + theOut->GetN()/2]);
      }
      else 
      {
        tempOut->SetPoint(i, theOut->GetX()[i], theOut->GetY()[i - theOut->GetN()/2]);
      }
    }
    delete theOut;
    AnalysisWaveform* tawf = new AnalysisWaveform(tempOut->GetN(), tempOut->GetY(), tempOut->GetX()[1] - tempOut->GetX()[0], 0);
    //tawf->padEven(1);
    theOut = new TGraph(tawf->phase()->GetN(), tawf->phase()->GetX(), tawf->phase()->GetY());
    for(int i = 0; i < theOut->GetN(); i++)
    {
      theOut->GetY()[i] *= TMath::RadToDeg();
    }
  }

  return theOut;
}

TH2D* PolarityMachine::runPolaritySimulation(int N, int eventNumber, double cr_phi, double cr_theta, double cr_snr, int pol, int metric, const char* outhistoname)
{
  d->getEvent(eventNumber, true);
  //make sure templates are loaded 
  loadTemplates(d->header()->realTime);
  FilterStrategy strat;
  TGraph* gCorr = 0;
  std::vector<double> corrs(numCRTemplates);

  FilteredAnitaEvent fae(d->useful(), &strat, d->gps(), d->header());
  wfcomb.combine(cr_phi, cr_theta, &fae, AnitaPol::AnitaPol_t(pol), 0, -25, 125);
  TGraph* g = new TGraph(wfcomb.getCoherent()->even()->GetN(), wfcomb.getCoherent()->even()->GetX(), wfcomb.getCoherent()->even()->GetY());
  g = FFTtools::normalizeWaveform(g);
  AnalysisWaveform* theCR = new AnalysisWaveform(g->GetN(), g->GetX(), g->GetY(), dT);
  delete g;
  TGraph* gg = new TGraph(theCR->even()->GetN(), theCR->even()->GetX(), theCR->even()->GetY());
  double t0 = gg->GetX()[0];
  for(int i = gg->GetN(); i < 2048; i++) gg->SetPoint(i, t0 + i*dT, 0); //pad to 2048
  
  for(int i = 0; i < numCRTemplates; i++)
  {
    g = new TGraph(simulatedCRs[i]->even()->GetN(), simulatedCRs[i]->even()->GetX(), simulatedCRs[i]->even()->GetY());
    gCorr = FFTtools::getCorrelationGraph(g, gg);
    double normalization = 1./(g->GetRMS(2) * gg->GetRMS(2) * g->GetN()/TMath::Power(2, int(TMath::Log2(g->GetN()))));
    double maxCorr = normalization * TMath::MaxElement(gCorr->GetN(), gCorr->GetY());
    double minCorr = normalization * TMath::MinElement(gCorr->GetN(), gCorr->GetY());
    corrs.at(i) = TMath::Max(maxCorr, TMath::Abs(minCorr));
    delete gCorr;
    delete g;
  }
  delete gg;
  std::vector<int> which_templates;
  double max_val = 0;
  int max_ind = -1;
  for(int j = 0; j < 3; j++)
  {
    max_val = -1;
    max_ind = -1;
    for(int i = 0; i < corrs.size(); i++)
    {
      if(max_val < corrs.at(i)) 
      {
        max_val = corrs.at(i);
        max_ind = i;
      }
    }
    which_templates.push_back(max_ind);
    corrs.erase(corrs.begin() + max_ind);
  }

  TH2D* hOut = generatePolarityMeasurements(N, eventNumber, cr_phi, cr_theta, cr_snr, pol, metric, outhistoname, which_templates);
  return hOut;
}

void PolarityMachine::makeSameSize(AnalysisWaveform* wf1, AnalysisWaveform* wf2)
{
  int n_target = (wf1->Neven() > wf1->Neven()) ? wf1->Neven() : wf2->Neven();
  if(wf1->Neven() != n_target)
  {
    TGraph g(n_target+2);
    for(int i = 0; i < g.GetN(); i++)
    {
      g.SetPoint(i, wf1->even()->GetX()[i], wf1->even()->GetY()[i]);
    }
    //delete wf1;
    wf1 = new AnalysisWaveform(g.GetN(), g.GetX(), g.GetY(), wf2->deltaT());
  }
  if(wf2->Neven() != n_target)
  {
    TGraph g(n_target+2);
    for(int i = 0; i < g.GetN(); i++)
    {
      g.SetPoint(i, wf2->even()->GetX()[i], wf2->even()->GetY()[i]);
    }
    //delete wf2;
    wf2 = new AnalysisWaveform(g.GetN(), g.GetX(), g.GetY(), wf1->deltaT());
  }
}


