#include "WaveformCombiner.h" 
#include "AntennaPositions.h" 
#include <vector> 
#include "BasicFilters.h" 
#include "FilteredAnitaEvent.h"
#include "AnitaVersion.h" 
#include "FFTtools.h"
#include "ResponseManager.h" 
#include "SystemResponse.h"
#include "DeltaT.h"
#include <dirent.h>
#include <sys/types.h>
#include <stdio.h>

UCorrelator::WaveformCombiner::WaveformCombiner(int nantennas, int npad, bool useUnfiltered, bool deconvolve, const AnitaResponse::ResponseManager * response, bool alfa_hack)
  : coherent(260), deconvolved(260), alfa_hack(alfa_hack)
{
  setNAntennas(nantennas); 
  setNPad(npad); 
  setDeconvolve(deconvolve); 
  setUseUnfiltered(useUnfiltered); 
  setGroupDelayFlag(true); 
  use_raw= useUnfiltered; 
  setResponseManager(response); 
  setBottomFirst(false);
  setDelayToCenter(false);
	extra_filters = 0;
	extra_filters_deconvolved = 0;
}



 

UCorrelator::WaveformCombiner::~WaveformCombiner()
{
}

const AnalysisWaveform * UCorrelator::WaveformCombiner::getDeconvolved() const 
{

  if (!do_deconvolution)
  {
//    fprintf(stderr,"WARNING! Deconvolution has not been enabled!\n"); 
    return 0; 
  }

  return &deconvolved; 
}




static void scaleGraph(TGraph * g, double C)
{

  int max = g->GetN(); 

  for (int i = 0; i < max; i++) 
  {
    g->GetY()[i] *= C; 
  }
}

static void addToGraph(TGraph * g, const TGraph * h)
{
  int max = TMath::Min(g->GetN(), h->GetN()); 
  for (int i = 0; i < max; i++)
  {
    g->GetY()[i] += h->GetY()[i]; 
  }
}


const AnalysisWaveform * UCorrelator::WaveformCombiner::wf (const FilteredAnitaEvent * event, int ant, AnitaPol::AnitaPol_t pol) 
{
  return use_raw ? event->getRawGraph(ant, pol) : event->getFilteredGraph(ant,pol); 


}


static SimplePassBandFilter alfa_filter(0,0.6);  

void UCorrelator::WaveformCombiner::combine(double phi, double theta, const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t pol, ULong64_t disallowed)
{

  std::vector<AnalysisWaveform> padded(nant);
  const int numDeconv = do_deconvolution ? nant : 0;
  std::vector<AnalysisWaveform> deconv(numDeconv);

  const UCorrelator::AntennaPositions * antpos = UCorrelator::AntennaPositions::instance(); 
  nant = antpos->getClosestAntennas(phi, nant, &antennas[0], disallowed); 
  double delays[nant]; 

  if (bottom_first) {
    if (antennas[0] < 32) {
      int toMove = antennas[0];
      if (antennas[1] > 32) {
       antennas[0] = antennas[1];
       antennas[1] = toMove; }
      else {
       antennas[0] = antennas[2];
       antennas[2] = toMove; }
    }
  }

  for (int i = 0; i < nant; i++) 
  {
    //ensure transform already calculated so we don't have to repeat when deconvolving
    (void) wf(event,antennas[i],pol)->freq(); 
		int ipol = (pol == AnitaPol::kVertical) ? 1 : 0;

    padded[i].~AnalysisWaveform(); 
    new (&padded[i]) AnalysisWaveform(*wf(event,antennas[i],pol));
		if(extra_filters)
		{
			for(int j = 0; j < extra_filters->nOperations(); j++)
			{
				FilterOperation * filterOp = (FilterOperation*) extra_filters->getOperation(j);
				filterOp->processOne(&padded[i], event->getHeader(), antennas[i], ipol);
				//delete filterOp;
			}
		}

    if (i == 0)
    {
      const TGraphAligned *p = wf(event,antennas[0],pol)->power(); 
      coherent_avg_spectrum.adopt(p); 
    }
    else
    {
      addToGraph(&coherent_avg_spectrum,wf(event,antennas[0],pol)->power()); 
    }

    //ALFA HACK IS HERE 
    if (alfa_hack && AnitaVersion::get() == 3 &&  pol == AnitaPol::kHorizontal && (antennas[i] == 4 || antennas[i] == 12))
    {
      alfa_filter.processOne(&padded[i]); 
    }

    padded[i].padFreq(npad);

    if (do_deconvolution)
    {
			deconv[i].~AnalysisWaveform(); 
      new (&deconv[i]) AnalysisWaveform(*wf(event,antennas[i],pol));
			if(extra_filters_deconvolved)
			{
				for(int j = 0; j < extra_filters_deconvolved->nOperations(); j++)
				{
					FilterOperation * filterOp = (FilterOperation*) extra_filters_deconvolved->getOperation(j);
					filterOp->processOne(&deconv[i], event->getHeader(), antennas[i], ipol);
					//delete filterOp;
				}
			}
      responses->response(pol,antennas[i])->deconvolveInPlace(&deconv[i], responses->getDeconvolutionMethod(), theta ); //TODO add angle  
      if (i == 0)
      {
        const TGraphAligned *p =deconv[i].power(); 
        deconvolved_avg_spectrum.adopt(p); 
      }
      else
      {
        addToGraph(&deconvolved_avg_spectrum, deconv[i].power()); 
      }
      deconv[i].padFreq(npad); 
    }

    if (delay_to_center) delays[i] = getDeltaTtoCenter(antennas[i], phi, theta, pol, enable_group_delay);
    else delays[i] = i == 0 ? 0 : getDeltaT(antennas[i], antennas[0], phi, theta, pol, enable_group_delay); 

  }

  combineWaveforms(nant, &padded[0], delays,0, &coherent); 
  scaleGraph(&coherent_avg_spectrum, 1./nant); 
  if (do_deconvolution)
  {
    combineWaveforms(nant, &deconv[0], delays,0, &deconvolved); 
    scaleGraph(&deconvolved_avg_spectrum, 1./nant); 
  }
}


AnalysisWaveform * UCorrelator::WaveformCombiner::combineWaveforms(int nwf, const AnalysisWaveform * wf, const double * delays, const double * scales, AnalysisWaveform * out)
{
  // we want to make the waveform as big as the entire span of all the waveforms to combine 

  double dt = wf[0].deltaT(); 


  double min = wf[0].even()->GetX()[0] + delays[0]; 
  double max = wf[0].even()->GetX()[wf[0].Neven()-1] + delays[0]; 

  for (int i = 1; i < nwf; i++) 
  {
    min = TMath::Min(min, wf[i].even()->GetX()[0] + delays[i]); 
    max = TMath::Max(max, wf[i].even()->GetX()[wf[i].Neven()-1] + delays[i]); 
  }

  int N = ceil((max-min)/dt); 

  if (!out) out = new AnalysisWaveform(N); 
  TGraph * gupdate = out->updateEven(); 
  gupdate->Set(N); 

  for (int i = 0; i < N; i++) 
  {
    gupdate->GetX()[i] = min + i * dt; 

    double val = 0; 
    for (int j = 0; j < nwf; j++) 
    {
      double scale = scales ? scales[j] : 1; 
      val += wf[j].evalEven(gupdate->GetX()[i] + delays[j])/nwf*scale; 
    }

//    printf("%f %f\n", gupdate->GetX()[i], val); 
    gupdate->GetY()[i] = val; 
  }

  return out; 
}

void UCorrelator::WaveformCombiner::setExtraFilters(FilterStrategy* extra)
{
	if(extra_filters) delete extra_filters;
	extra_filters = extra;
}

void UCorrelator::WaveformCombiner::setExtraFiltersDeconvolved(FilterStrategy* extra)
{
	if(extra_filters_deconvolved) delete extra_filters_deconvolved;
	extra_filters_deconvolved = extra;
}

