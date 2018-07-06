#include "FFTtools.h" 
#include "FilterStrategy.h"
#include "BH13Filter.h" 
#include "BasicFilters.h" 
#include "ResponseManager.h" 
#include "SystemResponse.h" 
#include "TGraph.h" 
#include "AnalysisWaveform.h" 
#include "RawAnitaHeader.h"
#include "TRandom.h"
#include "TFile.h"
#include "Adu5Pat.h"
#include "DigitalFilter.h"
#include "RawAnitaHeader.h"
#include <math.h>// for isnan
#include "AntennaPositions.h"
#include <ctype.h>
#include <stdio.h>

UCorrelator::BH13Filter::BH13Filter()
{
	TFile f(Form("%s/share/AnitaAnalysisFramework/responses/BH13TransferFn.root", getenv("ANITA_UTIL_INSTALL_DIR")));
	gPhase = (TGraph*) f.Get("fixPhase");
	gMag = (TGraph*) f.Get("fixAmp");
  f.Close();
}

UCorrelator::BH13Filter::~BH13Filter()
{
  delete gPhase;
  delete gMag;
}

void UCorrelator::BH13Filter::process(FilteredAnitaEvent * ev)
{
	for(int i = 0; i < 48; i++)
	{
		processOne(getWf(ev, i, AnitaPol::kHorizontal), ev->getHeader(), i, AnitaPol::kHorizontal);
		processOne(getWf(ev, i, AnitaPol::kVertical), ev->getHeader(), i, AnitaPol::kVertical);
	}
}

void UCorrelator::BH13Filter::processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol)
{
	if(whichAnt != 44 || whichPol != 0) return;
	AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal; 
  int ant = 44; 
	int old_size = awf->Neven();
	int nf = awf->Nfreq();
	double df = awf->deltaF();
	for( int i =0; i < nf; i++)
	{
		double f =i*df*1e9;
		double phase = awf->updateFreq()[i].getPhase()+gPhase->Eval(f);
		double mag = awf->updateFreq()[i].getAbs()*gMag->Eval(f);
		if( f>=.1e9 && f<=1.3e9) awf->updateFreq()[i].setMagPhase(mag, phase);
	}
	awf->updateEven()->Set(old_size);
}

UCorrelator::AntiBH13Filter::AntiBH13Filter()
{
	TFile f(Form("%s/share/AnitaAnalysisFramework/responses/BH13TransferFn.root", getenv("ANITA_UTIL_INSTALL_DIR")));
	gPhase = (TGraph*) f.Get("fixPhase");
	gMag = (TGraph*) f.Get("fixAmp");
  f.Close();
}

UCorrelator::AntiBH13Filter::~AntiBH13Filter()
{
  delete gPhase;
  delete gMag;
}

void UCorrelator::AntiBH13Filter::process(FilteredAnitaEvent * ev)
{
	for(int i = 0; i < 48; i++)
	{
		processOne(getWf(ev, i, AnitaPol::kHorizontal), ev->getHeader(), i, AnitaPol::kHorizontal);
		processOne(getWf(ev, i, AnitaPol::kVertical), ev->getHeader(), i, AnitaPol::kVertical);
	}
}

void UCorrelator::AntiBH13Filter::processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol)
{
	if(whichAnt != 44 || whichPol != 0) return;
	AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal; 
  int ant = 44; 
	int old_size = awf->Neven();
	int nf = awf->Nfreq();
	double df = awf->deltaF();
	for( int i =0; i < nf; i++)
	{
		double f =i*df*1e9;
		double phase = awf->updateFreq()[i].getPhase()-gPhase->Eval(f);
		double mag = awf->updateFreq()[i].getAbs()*(1./gMag->Eval(f));
		if( f>=.1e9 && f<=1.3e9) awf->updateFreq()[i].setMagPhase(mag, phase);
	}
	awf->updateEven()->Set(old_size);
}

void UCorrelator::timePadFilter::process(FilteredAnitaEvent* ev)
{
	for (int i = 0; i < 2*NUM_SEAVEYS; i++)
	{
		AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(i%2);
		int ant = i/2;
		getWf(ev, ant, pol)->padEven(fSamples);
	}
}

UCorrelator::A3toA4ConversionFilter::A3toA4ConversionFilter(char config)
{
  config = toupper(config);
  TFile f(Form("%s/share/AnitaAnalysisFramework/tuffModels/config%c.root", getenv("ANITA_UTIL_INSTALL_DIR"), config));
  gReal = (TGraph*) f.Get("gReal");
  gImag = (TGraph*) f.Get("gImag");
  phaseShift = TMath::ATan2(gImag->Eval(0), gReal->Eval(0));
  f.Close(); 
}

UCorrelator::A3toA4ConversionFilter::~A3toA4ConversionFilter()
{
  delete gReal;
  delete gImag;
}

void UCorrelator::A3toA4ConversionFilter::process(FilteredAnitaEvent* ev)
{
	for(int i = 0; i < 48; i++)
	{
		processOne(getWf(ev, i, AnitaPol::kHorizontal));
		processOne(getWf(ev, i, AnitaPol::kVertical));
	}
}

void UCorrelator::A3toA4ConversionFilter::processOne(AnalysisWaveform* awf, const RawAnitaHeader* header, int whichAnt, int whichPol)
{
	int old_size = awf->Neven();
	int nf = awf->Nfreq();
	double df = awf->deltaF();
  FFTWComplex* freq = awf->updateFreq();

	for( int i =0; i < nf; i++)
	{
    if(df*i > 1.5) continue;
    double reTemp = gReal->Eval(df*i);
    double imTemp = gImag->Eval(df*i);
    FFTWComplex temp(reTemp, imTemp);
    temp.setMagPhase(temp.getAbs(), temp.getPhase() - phaseShift);
    freq[i] *= temp;
	}
	awf->updateEven()->Set(old_size);
}






