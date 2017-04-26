#include "FFTtools.h" 
#include "FilterStrategy.h"
#include "BH13Filter.h" 
#include "SpectrumAverage.h" 
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
#include <math.h>// for isnan
#include "AntennaPositions.h"

void UCorrelator::BH13Filter::process(FilteredAnitaEvent * ev) 
{
	AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal; 
  int ant = 44; 
  AnalysisWaveform* awf = getWf(ev,ant,pol); 
	int old_size = awf->Neven();
	int nf = awf->Nfreq();
	double df = awf->deltaF();
	TFile f(Form("%s/share/UCorrelator/responses/BH13TransferFn.root", getenv("ANITA_UTIL_INSTALL_DIR")));
	TGraph* gPhase = (TGraph*) f.Get("fixPhase");
	TGraph* gMag = (TGraph*) f.Get("fixAmp");;
	for( int i =0; i < nf; i++)
	{
		double f =i*df*1e9;
		double phase = awf->updateFreq()[i].getPhase()+gPhase->Eval(f);
		double mag = awf->updateFreq()[i].getAbs()*gMag->Eval(f);
		if( f>=.1e9 && f<=1.3e9) awf->updateFreq()[i].setMagPhase(mag, phase);
	}
	awf->updateEven()->Set(old_size);
	delete gPhase;
	delete gMag;
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










