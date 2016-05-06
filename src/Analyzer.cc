#include "Analyzer.h" 
#include "AnalysisConfig.h" 
#include "AnalysisWaveform.h"
#include "FilteredAnitaEvent.h"
#include "Util.h"
#include "DigitalFilter.h" 
#include "FFTtools.h" 
#include "UsefulAdu5Pat.h"
#include "PeakFinder.h"
#include "Flags.h"


static UCorrelator::AnalysisConfig defaultConfig; 
static int instance_counter = 0; 


UCorrelator::Analyzer::Analyzer(const AnalysisConfig * conf, bool interactive) 
  : cfg(conf ? conf: &defaultConfig),
    corr(cfg->correlator_nphi,0,360,  cfg->correlator_ntheta, -cfg->correlator_theta_lowest, cfg->correlator_theta_highest) , 
    wfcomb(cfg->combine_nantennas, cfg->combine_npad), 
    zoomed(TString::Format("zoomed_%d", instance_counter), "Zoomed!", cfg->zoomed_nphi, 0 ,1, cfg->zoomed_ntheta, 0, 1),
   interactive(interactive)  
{

  corr.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb.setGroupDelayFlag(cfg->enable_group_delay); 
  instance_counter++; 
  power_filter = 0; //TODO initialize this 

  if (interactive) 
  {
    correlation_maps[0] = new TH2D; 
    correlation_maps[1] = new TH2D; 

    for (int i = 0; i < cfg->nmaxima; i++)
    {
      zoomed_correlation_maps[0].push_back(new TH2D); 
      zoomed_correlation_maps[1].push_back(new TH2D); 
      coherent[0].push_back(new AnalysisWaveform); 
      coherent[1].push_back(new AnalysisWaveform); 
      deconvolved[0].push_back(new AnalysisWaveform); 
      deconvolved[1].push_back(new AnalysisWaveform); 

    }

  }
}




void UCorrelator::Analyzer::analyze(const FilteredAnitaEvent * event, AnitaEventSummary * summary) 
{
//  summary->~AnitaEventSummary(); 
  summary = new (summary) AnitaEventSummary(event->getHeader(), (UsefulAdu5Pat*) event->getGPS()); 

  //check for saturation
  uint64_t saturated[2] = {0,0}; 
  UCorrelator::flags::checkSaturation(event->getUsefulAnitaEvent(), 
                                      &saturated[AnitaPol::kHorizontal], 
                                      &saturated[AnitaPol::kVertical], 
                                      cfg->saturation_threshold); 

  //we need a UsefulAdu5Pat for this event
  
  UsefulAdu5Pat * pat =  (UsefulAdu5Pat*) event->getGPS();  //unconstifying it .. hopefully that won't cause problems

  // These will be needed for peak finding

  // loop over wanted polarizations 


  for (int pol = cfg->start_pol; pol <= cfg->end_pol; pol++) 
  {
    // tell the correlator not to use saturated events and make the correlation map
    corr.setDisallowedAntennas(saturated[pol]); 
    corr.compute(event, AnitaPol::AnitaPol_t(pol)); 

    //compute RMS of correlation map 
    maprms = corr.getHist()->GetRMS(3); 

     if (interactive) 
     {
       correlation_maps[pol]->~TH2D(); 
       correlation_maps[pol] = new (correlation_maps[pol]) TH2D(* (TH2D*) corr.getHist()); 
     }

    // Find the isolated peaks in the image 
    peakfinder::RoughMaximum maxima[cfg->nmaxima]; 
    int npeaks = UCorrelator::peakfinder::findIsolatedMaxima((const TH2D*) corr.getHist(), cfg->peak_isolation_requirement, cfg->nmaxima, maxima, cfg->use_bin_center); 
//    printf("npeaks: %d\n", npeaks); 
    summary->nPeaks[pol] = npeaks; 

    // Loop over found peaks 
    for (int i = 0; i < npeaks; i++) 
    {
      // zoom in on the values 
      printf("rough phi:%f, rough theta: %f\n", maxima[i].x, -maxima[i].y); 
      fillPointingInfo(maxima[i].x, maxima[i].y, &summary->peak[pol][i], pat); 
      if (interactive)
      {
        zoomed_correlation_maps[pol][i]->~TH2D(); 
        zoomed_correlation_maps[pol][i] = new (zoomed_correlation_maps[pol][i]) TH2D(zoomed); 
      }
      printf("phi:%f, theta:%f\n", summary->peak[pol][i].phi, summary->peak[pol][i].theta); 
      //now make the combined waveforms 
      wfcomb.combine(summary->peak[pol][i].phi, -summary->peak[pol][i].theta, event, (AnitaPol::AnitaPol_t) pol, saturated[pol]); 
      fillWaveformInfo(wfcomb.getCoherent(), &summary->coherent[pol][i]); 
      fillWaveformInfo(wfcomb.getDeconvolved(), &summary->deconvolved[pol][i]); 

      if (interactive)
      {
        coherent[pol][i]->~AnalysisWaveform(); 
        deconvolved[pol][i]->~AnalysisWaveform(); 
        coherent[pol][i] = new (coherent[pol][i]) AnalysisWaveform(*wfcomb.getCoherent()); 
        deconvolved[pol][i] = new (deconvolved[pol][i]) AnalysisWaveform(*wfcomb.getDeconvolved()); 
      }
    }
  }

  fillFlags(event, &summary->flags, pat); 
}

void UCorrelator::Analyzer::fillPointingInfo(double rough_phi, double rough_theta, AnitaEventSummary::PointingHypothesis * point, UsefulAdu5Pat * pat)
{
      corr.computeZoomed(rough_phi, rough_theta, cfg->zoomed_nphi, cfg->zoomed_dphi,  cfg->zoomed_ntheta, cfg->zoomed_dtheta, cfg->zoomed_nant, &zoomed); 

      //get pointer to the pointing hypothesis we are about to fill 

      // This will fill in phi, theta, value, var_theta, var_phi and covar 
      
      peakfinder::FineMaximum max; 
      switch (cfg->fine_peak_finding_option)
      {
        case AnalysisConfig::FinePeakFindingAbby: 
          UCorrelator::peakfinder::doInterpolationPeakFindingAbby(&zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingBicubic: 
          UCorrelator::peakfinder::doInterpolationPeakFindingBicubic(&zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingQuadraticFit16: 
          UCorrelator::peakfinder::doPeakFindingQuadratic16(&zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingQuadraticFit25: 
          UCorrelator::peakfinder::doPeakFindingQuadratic25(&zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingGaussianFit: 
          UCorrelator::peakfinder::doPeakFindingGaussian(&zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingQuadraticFit9: 
        default: 
          UCorrelator::peakfinder::doPeakFindingQuadratic9(&zoomed, &max); 
          break; 
      }; 
      
      max.copyToPointingHypothesis(point); 

      //snr is ratio of point value to map rms
      point->snr = point->value / maprms; 

      // will do this later
      point->hwAngle = -9999; // TODO 

      //Compute intersection with continent, or set values to -9999 if no intersection
      if (!pat->getSourceLonAndLatAtAlt(point->phi, point->theta, point->latitude, point->longitude, point->altitude)) 
      {
        point->latitude = -9999; 
        point->longitude = -9999;  
        point->altitude = -9999; 
      }
}


void UCorrelator::Analyzer::fillWaveformInfo(const AnalysisWaveform * wf, AnitaEventSummary::WaveformInfo * info)
{
  if (!wf)
  {
    memset(info, 0, sizeof(AnitaEventSummary::WaveformInfo)); 
    return; 
  }
  const TGraph * even = wf->even(); 
  info->peakVal = FFTtools::getPeakVal((TGraph*) even); 
  info->peakHilbert = FFTtools::getPeakVal((TGraph*) wf->hilbertEnvelope()); 
  info->numAntennasInCoherent = cfg->combine_nantennas; 

  double dt = wf->deltaT(); 
  double t0 = even->GetX()[0]; 

  int i0 = TMath::Max(0.,floor((cfg->noise_estimate_t0 - t0)/dt)); 
  int i1 = TMath::Min(even->GetN()-1.,ceil((cfg->noise_estimate_t1 - t0)/dt)); 
  int n = i1 - i0 + 1; 
//  printf("%d-%d -> %d \n", i0, i1, n); 

  double rms = TMath::RMS(n, even->GetY() + i0); 
  
  info->snr = info->peakVal / rms; 

  TGraph power(*wf->powerdB()); 

  if (power_filter)
  {
    power_filter->filterGraph(&power); 
  }

  int peak_index; 
  double peak_val =  FFTtools::getPeakVal(&power,&peak_index); 
  int imax = power.GetN(); 
  int imin = 0; 


  for (int i = peak_index+1; i <power.GetN(); i++)
  {
    if (power.GetY()[i] < peak_val - cfg->bw_ndb)
    {
      imax = i; 
      break; 
    }
  }

  for (int i = peak_index+1; i >= 0; i--)
  {
    if (power.GetY()[i] < peak_val - cfg->bw_ndb)
    {
      imin = i; 
      break; 
    }
  }

  info->bandwidth = (imax - imin) * (power.GetX()[1]-power.GetX()[0]); 


}


UCorrelator::Analyzer::~Analyzer()
{
  if (interactive)
  {
    delete correlation_maps[0];
    delete correlation_maps[1]; 

    for (int i = 0; i < cfg->nmaxima; i++)
    {
      delete zoomed_correlation_maps[0][i];
      delete zoomed_correlation_maps[1][i];
      delete coherent[0][i];
      delete coherent[1][i];
      delete deconvolved[0][i];
      delete deconvolved[1][i];
    }
  }
}

void UCorrelator::Analyzer::fillFlags(const FilteredAnitaEvent * fae, AnitaEventSummary::EventFlags * flags, UsefulAdu5Pat * pat) 
{

  //TODO 
  flags->isVarner = false; 
  flags->isVarner2 = false; 
  flags->nadirFlag = true; 
  flags->strongCWFlag = false; 


  if ( isLDBHPol(pat, fae->getHeader(), cfg) || isLDBVPol (pat, fae->getHeader(), cfg))
  {
    flags->pulser = AnitaEventSummary::EventFlags::LDB; 
  }
  else if ( isWAISHPol(pat, fae->getHeader(), cfg) || isWAISVPol (pat, fae->getHeader(), cfg))
  {
    flags->pulser = AnitaEventSummary::EventFlags::WAIS; 
  }
  else
  {
    flags->pulser = AnitaEventSummary::EventFlags::NONE; 
  }

  flags->isGood = !flags->isVarner && !flags->isVarner2 && !flags->strongCWFlag; 

}




