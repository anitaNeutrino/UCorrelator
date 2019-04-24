#ifndef POLARITY_MACHINE_H
#define POLARITY_MACHINE_H

#include "TObject.h"
#include "FFTtools.h"
#include "TFile.h"
#include "AnalysisWaveform.h"
#include "AnitaVersion.h"
#include "AnitaDataset.h"
#include "TRandom3.h"
#include "RawAnitaHeader.h"
#include "WaveformCombiner.h"
#include "FilteredAnitaEvent.h"
#include "TH2D.h"


/**
  This is a class meant to explore polarity measurements.
  For now it's just for EAS-like events, but someday will maybe have neutrino-like events, as well
  Andrew Ludwig <abl@uchicago.edu>
  */

namespace AnitaResponse
{
  class DeconvolutionMethod; 
}


class PolarityMachine
{
  public:
    PolarityMachine(int padX=15);

    ~PolarityMachine();
    void zeroInternals();

    void loadTemplates(unsigned int eventTime, int version = AnitaVersion::get());

    /* N is number of times you want to do this
     * eventNumber is the evetn you want to generate the noise using the min bias triggers around
     * cr_phi, cr_theta, cr_snr are the properties of the event you are modeling (payload coordinates)
     * pol is 0 or 1 for H or V
     * metric is how you want to classify polarity
     * metric is also used in the testPolarity function - corresponding numbers are printed by that fn
     * outfilename is the file where you want to save histograms to.
     * whichTemplates is which templates you want to use to check this (default is 0-4) */

    TH2D* generatePolarityMeasurements(int N, int eventNumber, double cr_phi = 0, double cr_theta = 0, double cr_snr = 3, int pol = 0, int metric = 0,  const char * outhistoname = "h0", std::vector<int> whichTemplates={0,1,2,3,4});

    static const int numCRTemplates = 10;
    AnalysisWaveform* simulatedCRs[numCRTemplates];
    double simulatedCRPeaks[numCRTemplates];

    AnalysisWaveform* theImpTemplate;
    AnalysisWaveform* CRtemplate;
    AnalysisWaveform* CRtemplate_deconv;
    AnalysisWaveform* CRtemplate_deconv_windowed;

    void setDeconvolutionMethod(AnitaResponse::DeconvolutionMethod* opt) { deconv = opt ;}
    void setWindowSizeNs(double opt) { windowSizeNs = opt ;}
    void setPadFactor(int opt) { padFactor = opt ;}
    void setTestCoherent(bool opt) { test_coherent = opt ;}
    void setOffset(int opt) { offset = opt ;}

    /* Metrics are as follows:
     *    0 is largest peak in the cross correlation template being negative or positive (binary)
     *    1 is peak to sidelobe ratio
     *    2 same as 0 but w/ a window about the deconvolved wf
     *    3 same as 1 but w/ a window about the deconvolved wf
     *    4 is Andres' metric of (max+min)/(max-min) of wf (no correlation)
     *    5 is the Peter metric (windows differently and does a fixed point correlation about a fiducial)
     *    6 is the fourier phase metric
     * Will add more as i think of them */

    double testPolarity(int metric, AnalysisWaveform* wf, bool deconvolved);
   
    
    /* wf is the waveform you want to test.
     * Windowing types are as follows:
     *    0 windows using Peter's fiducial (zero crossing nearest max stokes I)
     *    1 picks its fiducial at the mean of the CDF of the hilbert peak
     *    2 picks its fiducial at the 90% max of the CDF of the hilbert peak
     *    3 picks its fiducial at the max of the hilbert peak
     * Will add more as i think of them 
     * returns the average value of the phase */
    double fourierPhasePolarity(int window_type, AnalysisWaveform* wf);

    TH2D* runPolaritySimulation(int N , int eventNumber, double cr_phi = 0, double cr_theta = 0, double cr_snr = 3, int pol = 0, int metric = 0, const char * outhistoname = "h0");

    /* draws a nearby min bias event */
    AnalysisWaveform* makeNoiseWaveformFromMinBias(int eventNumber, double phi, double theta, int pol, int current_N);

    /* adds the nearby min bias event to the template waveform with some sort of snr scaling */
    AnalysisWaveform* generateNoisyWaveform(AnalysisWaveform* templateWf, AnalysisWaveform* noiseWf, double snr, TRandom3* tr);

    /* windows the waveform about the maximum stokes I (tukey window, not a trimmed wf) */
    TGraph* windowWaveform(AnalysisWaveform* wf, double window_size_ns = 5);

    /* windows the waveform about the zero crossing with the maximum stokes I */
    TGraph* peterWindow(AnalysisWaveform* wf, double window_size_ns = 5, bool roll = false);

    /* windows the waveform using something about the hilbert peak
     * wf is the input waveform
     * window_size_ns is the total size of the output you want
     * roll means do you want to roll the fiducial to the front (yes for fourier phase method) 
     * use_cdf is a bool to use the cdf of the hilbert peak to pick the fiducial 
     * how_to_use_cdf 0 = mean of cdf, 1 = 90% max of cdf
     */
    TGraph* hilbertWindow(AnalysisWaveform* wf, double window_size_ns = 5, bool roll = false, bool use_cdf = false, int how_to_use_cdf = 0);

    /* makes two analysis waveforms the same size (useful for correlation) */
    void makeSameSize(AnalysisWaveform* wf1, AnalysisWaveform* wf2);

    /* these functions pull the two waveforms that went into the most recently formed correlation graph and the correlation graph itself, mostly for validation purposes */
    TGraph* getInputWaveform(){return inputWaveform ;}
    TGraph* getCRWaveform(){return crWaveform ;}
    TGraph* getCorrelationGraph(){return correlationGraph ;}

    std::string getNotchStr(){return fNotchStr ;}

  private:
    AnitaDataset* d;
    int length;
    int lengthFFT;
    int padFactor;
    double dT;
    double dF;
    double windowSizeNs;
    bool test_coherent;
    int forwardMinBias;
    int previousMinBias;
    int offset;

    bool kTemplatesLoaded;
    std::string fNotchStr;
    std::vector<int> payloadTimes;
    std::vector<std::string> notchConfigs;

    AnitaResponse::DeconvolutionMethod * deconv;
    UCorrelator::WaveformCombiner wfcomb;

    void fillNotchConfigs();
    void getCRTemplates(int version = AnitaVersion::get());
    void getImpulseResponseTemplate(int version = AnitaVersion::get());

    TGraph* inputWaveform;
    TGraph* crWaveform;
    TGraph* correlationGraph;


    ClassDefNV(PolarityMachine, 1);

};

#endif
