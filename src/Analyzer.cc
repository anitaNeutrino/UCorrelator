#include "Analyzer.h" 
#include "AnalysisConfig.h" 
#include "TPaveText.h"
#include "AnalysisWaveform.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "FilteredAnitaEvent.h"
#include "Util.h"
#include "DigitalFilter.h" 
#include "FFTtools.h" 
#include "UsefulAdu5Pat.h"
#include "RawAnitaHeader.h" 
#include "PeakFinder.h"
#include "UCFlags.h"
#include "ShapeParameters.h" 
#include "SpectrumParameters.h" 


static UCorrelator::AnalysisConfig defaultConfig; 
static int instance_counter = 0; 

#ifndef DEG2RAD
#define DEG2RAD (TMath::Pi()/ 180) 
#endif

#ifndef RAD2DEG
#define RAD2DEG (180 / TMath::Pi()) 
#endif






UCorrelator::Analyzer::Analyzer(const AnalysisConfig * conf, bool interactive) 
  : cfg(conf ? conf: &defaultConfig),
    corr(cfg->correlator_nphi,0,360,  cfg->correlator_ntheta, -cfg->correlator_theta_lowest, cfg->correlator_theta_highest) , 
    wfcomb(cfg->combine_nantennas, cfg->combine_npad), 
    wfcomb_xpol(cfg->combine_nantennas, cfg->combine_npad), 
    zoomed(TString::Format("zoomed_%d", instance_counter), "Zoomed!", cfg->zoomed_nphi, 0 ,1, cfg->zoomed_ntheta, 0, 1),
   interactive(interactive)  
{

  corr.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb_xpol.setGroupDelayFlag(cfg->enable_group_delay); 
  instance_counter++; 
  power_filter = new FFTtools::GaussianFilter(2,3) ; //TODO make this configurable

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
      coherent_xpol[0].push_back(new AnalysisWaveform); 
      coherent_xpol[1].push_back(new AnalysisWaveform); 
      deconvolved_xpol[0].push_back(new AnalysisWaveform); 
      deconvolved_xpol[1].push_back(new AnalysisWaveform); 

      coherent_power[0].push_back(new TGraphAligned); 
      coherent_power[1].push_back(new TGraphAligned); 
      deconvolved_power[0].push_back(new TGraphAligned); 
      deconvolved_power[1].push_back(new TGraphAligned); 
      coherent_power_xpol[0].push_back(new TGraphAligned); 
      coherent_power_xpol[1].push_back(new TGraphAligned); 
      deconvolved_power_xpol[0].push_back(new TGraphAligned); 
      deconvolved_power_xpol[1].push_back(new TGraphAligned); 
    }
  }
}





void UCorrelator::Analyzer::analyze(const FilteredAnitaEvent * event, AnitaEventSummary * summary) 
{
  
  /* Initialize the summary */ 
  summary = new (summary) AnitaEventSummary(event->getHeader(), (UsefulAdu5Pat*) event->getGPS()); 

  //check for saturation
  uint64_t saturated[2] = {0,0}; 
  UCorrelator::flags::checkSaturation(event->getUsefulAnitaEvent(), 
                                      &saturated[AnitaPol::kHorizontal], 
                                      &saturated[AnitaPol::kVertical], 
                                      cfg->saturation_threshold); 

  //we need a UsefulAdu5Pat for this event
  UsefulAdu5Pat * pat =  (UsefulAdu5Pat*) event->getGPS();  //unconstifying it .. hopefully that won't cause problems
  

  // loop over wanted polarizations 
  for (int pol = cfg->start_pol; pol <= cfg->end_pol; pol++) 
  {

    UShort_t triggeredPhi = AnitaPol::AnitaPol_t(pol) == AnitaPol::kHorizontal ? event->getHeader()->l3TrigPatternH : event->getHeader()->l3TrigPattern; 
    UShort_t triggeredPhiL1 = AnitaPol::AnitaPol_t(pol) == AnitaPol::kHorizontal ? event->getHeader()->l1TrigMaskH : event->getHeader()->l1TrigMask;  
    UShort_t maskedPhi = AnitaPol::AnitaPol_t(pol) == AnitaPol::kHorizontal ? event->getHeader()->phiTrigMaskH : event->getHeader()->phiTrigMask; 

    TVector2 triggerAngle; 

    int ntriggered = __builtin_popcount(triggeredPhi); 
    int ntriggeredL1 = ntriggered ? 0 : __builtin_popcount(triggeredPhiL1); 

    UShort_t which_trigger = ntriggered ? triggeredPhi : triggeredPhiL1; 
    int which_ntriggered = ntriggered ? ntriggered : ntriggeredL1; 

    for (int i = 0; i < NUM_PHI; i++) 
    {
      if (which_trigger & (1 << i))
      {
         double ang = (i * 22.5) * TMath::Pi()/180;
         triggerAngle += TVector2(cos(ang), sin(ang)) / which_ntriggered; 
      }
    }

    double avgHwAngle = triggerAngle.Phi(); 

    // tell the correlator not to use saturated events and make the correlation map
    corr.setDisallowedAntennas(saturated[pol]); 
    corr.compute(event, AnitaPol::AnitaPol_t(pol)); 

    //compute RMS of correlation map 
    maprms = corr.getHist()->GetRMS(3); 

    // Find the isolated peaks in the image 
    peakfinder::RoughMaximum maxima[cfg->nmaxima]; 
    int npeaks = UCorrelator::peakfinder::findIsolatedMaxima((const TH2D*) corr.getHist(), cfg->peak_isolation_requirement, cfg->nmaxima, maxima, cfg->use_bin_center); 
//    printf("npeaks: %d\n", npeaks); 
    summary->nPeaks[pol] = npeaks; 

    rough_peaks[pol].clear(); 

    // Loop over found peaks 
    for (int i = 0; i < npeaks; i++) 
    {
      // zoom in on the values 
//      printf("rough phi:%f, rough theta: %f\n", maxima[i].x, -maxima[i].y); 

      rough_peaks[pol].push_back(std::pair<double,double>(maxima[i].x, maxima[i].y)); 

      fillPointingInfo(maxima[i].x, maxima[i].y, &summary->peak[pol][i], pat, avgHwAngle, triggeredPhi, maskedPhi); 


      //fill in separation 
      summary->peak[pol][i].phi_separation = 1000; 
      for (int j = 1; j < i; j++)
      {
        summary->peak[pol][i].phi_separation = TMath::Min(summary->peak[pol][i].phi_separation, fabs(FFTtools::wrap(summary->peak[pol][i].phi - summary->peak[pol][j].phi, 360, 0))); 
      }

//      printf("phi:%f, theta:%f\n", summary->peak[pol][i].phi, summary->peak[pol][i].theta); 


      //now make the combined waveforms 
      wfcomb.combine(summary->peak[pol][i].phi, -summary->peak[pol][i].theta, event, (AnitaPol::AnitaPol_t) pol, saturated[pol]); 
      wfcomb_xpol.combine(summary->peak[pol][i].phi, -summary->peak[pol][i].theta, event, (AnitaPol::AnitaPol_t) (1-pol), saturated[pol]); 
      fillWaveformInfo(wfcomb.getCoherent(), wfcomb_xpol.getCoherent(), wfcomb.getCoherentAvgSpectrum(), &summary->coherent[pol][i], (AnitaPol::AnitaPol_t) pol); 
      fillWaveformInfo(wfcomb.getDeconvolved(), wfcomb_xpol.getDeconvolved(), wfcomb.getDeconvolvedAvgSpectrum(), &summary->deconvolved[pol][i],  (AnitaPol::AnitaPol_t)pol); 


      if (interactive) //copy everythig
      {

        zoomed_correlation_maps[pol][i]->~TH2D(); 
        zoomed_correlation_maps[pol][i] = new (zoomed_correlation_maps[pol][i]) TH2D(zoomed); 

        coherent[pol][i]->~AnalysisWaveform(); 
        coherent[pol][i] = new (coherent[pol][i]) AnalysisWaveform(*wfcomb.getCoherent()); 

        coherent_power[pol][i]->~TGraphAligned(); 
        coherent_power[pol][i] = new (coherent_power[pol][i]) TGraphAligned(*wfcomb.getCoherentAvgSpectrum()); 
        coherent_power[pol][i]->dBize(); 


        if (wfcomb.getDeconvolved())
        {
          deconvolved[pol][i]->~AnalysisWaveform(); 
          deconvolved[pol][i] = new (deconvolved[pol][i]) AnalysisWaveform(*wfcomb.getDeconvolved()); 
          deconvolved[pol][i]->updateEven()->SetLineColor(2); 

          deconvolved_power[pol][i]->~TGraphAligned(); 
          deconvolved_power[pol][i] = new (deconvolved_power[pol][i]) TGraphAligned(*wfcomb.getDeconvolvedAvgSpectrum()); 
          deconvolved_power[pol][i]->dBize(); 
          deconvolved_power[pol][i]->SetLineColor(2); 
          interactive_deconvolved = true; 
        }
        else
        {
          interactive_deconvolved = false; 
        }


        coherent_xpol[pol][i]->~AnalysisWaveform(); 
        coherent_xpol[pol][i] = new (coherent_xpol[pol][i]) AnalysisWaveform(*wfcomb_xpol.getCoherent()); 
        coherent_xpol[pol][i]->updateEven()->SetLineColor(11); 
        coherent_xpol[pol][i]->updateEven()->SetLineStyle(3); 
        
        coherent_power_xpol[pol][i]->~TGraphAligned(); 
        coherent_power_xpol[pol][i] = new (coherent_power_xpol[pol][i]) TGraphAligned(*wfcomb_xpol.getCoherentAvgSpectrum()); 
        coherent_power_xpol[pol][i]->dBize(); 
        coherent_power_xpol[pol][i]->SetLineStyle(3); 
        coherent_power_xpol[pol][i]->SetLineColor(11); 


        if (wfcomb_xpol.getDeconvolved())
        {
          deconvolved_xpol[pol][i]->~AnalysisWaveform(); 
          deconvolved_xpol[pol][i] = new (deconvolved_xpol[pol][i]) AnalysisWaveform(*wfcomb_xpol.getDeconvolved()); 

          deconvolved_xpol[pol][i]->updateEven()->SetLineColor(45); 
          deconvolved_xpol[pol][i]->updateEven()->SetLineStyle(3); 
          deconvolved_power_xpol[pol][i]->~TGraphAligned(); 
          deconvolved_power_xpol[pol][i] = new (deconvolved_power_xpol[pol][i]) TGraphAligned(*wfcomb_xpol.getDeconvolvedAvgSpectrum()); 
          deconvolved_power_xpol[pol][i]->dBize(); 
          deconvolved_power_xpol[pol][i]->SetLineStyle(3); 
          deconvolved_power_xpol[pol][i]->SetLineColor(46); 
          interactive_xpol_deconvolved = true; 
        }
        else
        {
          interactive_xpol_deconvolved = false; 
        }
      }
    }
    if (interactive) 
    {
       correlation_maps[pol]->~TH2D(); 
       correlation_maps[pol] = new (correlation_maps[pol]) TH2D(* (TH2D*) corr.getHist()); 
    }
  }

  fillFlags(event, &summary->flags, pat); 

  if (interactive) last = *summary; 

}

void UCorrelator::Analyzer::fillPointingInfo(double rough_phi, double rough_theta, AnitaEventSummary::PointingHypothesis * point, 
                                             UsefulAdu5Pat * pat, double hwPeakAngle, UShort_t triggered_sectors, UShort_t masked_sectors)
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

      point->hwAngle = FFTtools::wrap(point->phi - hwPeakAngle,360,0); 
      int sector = fmod(point->phi + 11.25,360) / 22.5; 

      point->masked = masked_sectors & ( 1 << sector); 
      point->triggered = triggered_sectors & ( 1 << sector); 

      //Compute intersection with continent, or set values to -9999 if no intersection
      if (!pat->traceBackToContinent(point->phi * DEG2RAD, point->theta * DEG2RAD, &point->longitude, &point->latitude, &point->altitude, &point->theta_adjustment_needed)) 
      {
        point->latitude = -9999; 
        point->longitude = -9999;  
        point->altitude = -9999; 
        point->theta_adjustment_needed = -9999; 
      }
      else
      {
        point->theta_adjustment_needed *= RAD2DEG; 
      }
}


void UCorrelator::Analyzer::fillWaveformInfo(const AnalysisWaveform * wf, const AnalysisWaveform * xpol_wf, const TGraph* pwr, AnitaEventSummary::WaveformInfo * info, AnitaPol::AnitaPol_t pol)
{
  if (!wf)
  {
    memset(info, 0, sizeof(AnitaEventSummary::WaveformInfo)); 
    return; 
  }
  const TGraphAligned * even = wf->even(); 
  const TGraphAligned * xpol_even= xpol_wf->even(); 
  info->peakVal = FFTtools::getPeakVal((TGraph*) even); 
  info->xPolPeakVal = FFTtools::getPeakVal((TGraph*) xpol_even); 
  info->peakHilbert = FFTtools::getPeakVal((TGraph*) wf->hilbertEnvelope()); 
  info->xPolPeakHilbert = FFTtools::getPeakVal((TGraph*) xpol_wf->hilbertEnvelope()); 
  info->numAntennasInCoherent = cfg->combine_nantennas; 

  info->totalPower = even->getSumV2(); 
  info->totalPowerXpol = xpol_even->getSumV2(); 

  info->riseTime_10_90 = shape::getRiseTime((TGraph*) wf->hilbertEnvelope(), 0.1*info->peakVal, 0.9*info->peakVal); 
  info->riseTime_10_50 = shape::getRiseTime((TGraph*) wf->hilbertEnvelope(), 0.1*info->peakVal, 0.5*info->peakVal); 
  info->fallTime_90_10 = shape::getFallTime((TGraph*) wf->hilbertEnvelope(), 0.1*info->peakVal, 0.9*info->peakVal); 
  info->fallTime_50_10 = shape::getFallTime((TGraph*) wf->hilbertEnvelope(), 0.1*info->peakVal, 0.5*info->peakVal); 

  int ifirst, ilast; 
  info->width_50_50 = shape::getWidth((TGraph*) wf->hilbertEnvelope(), 0.5*info->peakVal, &ifirst, &ilast); 
  info->power_50_50 = even->getSumV2(ifirst, ilast); 



  info->width_10_10 = shape::getWidth((TGraph*) wf->hilbertEnvelope(), 0.1*info->peakVal, &ifirst, &ilast); 
  info->power_10_10 = even->getSumV2(ifirst, ilast); 




  if (pol == AnitaPol::kHorizontal)
  {
    FFTtools::stokesParameters(even->GetN(), 
                               even->GetY(), 
                               wf->hilbertTransform()->even()->GetY(), 
                               xpol_even->GetY(), 
                               xpol_wf->hilbertTransform()->even()->GetY(), 
                               &(info->I), &(info->Q), &(info->U), &(info->V)); 
  }
  else
  {
    FFTtools::stokesParameters(even->GetN(), 
                               xpol_even->GetY(), 
                               xpol_wf->hilbertTransform()->even()->GetY(), 
                               even->GetY(), 
                               wf->hilbertTransform()->even()->GetY(), 
                               &(info->I), &(info->Q), &(info->U), &(info->V)); 
 
  }

  double dt = wf->deltaT(); 
  double t0 = even->GetX()[0]; 

  int i0 = TMath::Max(0.,floor((cfg->noise_estimate_t0 - t0)/dt)); 
  int i1 = TMath::Min(even->GetN()-1.,ceil((cfg->noise_estimate_t1 - t0)/dt)); 
  int n = i1 - i0 + 1; 
//  printf("%d-%d -> %d \n", i0, i1, n); 

  double rms = TMath::RMS(n, even->GetY() + i0); 
  
  info->snr = info->peakVal / rms; 

  TGraphAligned power(*pwr); 
  power.dBize(); 

  if (power_filter)
  {
    power_filter->filterGraph(&power); 
  }

  spectrum::fillSpectrumParameters(&power, info, cfg); 

}


UCorrelator::Analyzer::~Analyzer()
{
  if (interactive)
  {
    delete correlation_maps[0];
    delete correlation_maps[1]; 


    for (int pol = 0; pol < 2; pol++)
    {
      for (int i = 0; i < cfg->nmaxima; i++)
      {
        delete zoomed_correlation_maps[pol][i];

        delete coherent[pol][i];
        delete deconvolved[pol][i];

        delete coherent_xpol[pol][i];
        delete deconvolved_xpol[pol][i];

        delete coherent_power[pol][i];
        delete deconvolved_power[pol][i];

        delete coherent_power_xpol[pol][i];
        delete deconvolved_power_xpol[pol][i];
      }
    }

    clearInteractiveMemory(); 
  }

  if (power_filter)
    delete power_filter; 
}

void UCorrelator::Analyzer::clearInteractiveMemory() const
{

  for (unsigned i = 0; i < delete_list.size(); i++) 
  {
    delete delete_list[i]; 
  }

  delete_list.clear(); 
}


static void setOnClickHandler(TPad * pad) 
{

}





void UCorrelator::Analyzer::drawSummary(TPad * ch, TPad * cv) const
{
  TPad * pads[2] = {ch,cv}; 

  clearInteractiveMemory(); 

  gStyle->SetOptStat(0); 
  for (int ipol = cfg->start_pol; ipol <= cfg->end_pol; ipol++)
  {
    if (!pads[ipol])
    {
      pads[ipol] = new TCanvas(ipol == 0 ? "analyzer_ch" : "analyzer_cv", ipol == 0 ? "hpol" : "vpol",1920,500); 
    }

    pads[ipol]->Clear(); 
    pads[ipol]->Divide(2,1); 

    pads[ipol]->cd(1)->Divide(1,2); 

    pads[ipol]->cd(1)->cd(1); 
    correlation_maps[ipol]->SetTitle(ipol == 0 ? "HPol map" : "VPol map" ); 
    correlation_maps[ipol]->Draw("colz"); 

    for (int i = 0; i < last.nPeaks[ipol]; i++) 
    {
      pads[ipol]->cd(1)->cd(1); 
      TMarker * m = new TMarker(rough_peaks[ipol][i].first, rough_peaks[ipol][i].second, 3); 
      m->SetMarkerSize(last.nPeaks[ipol] -i); 
      m->Draw(); 
      delete_list.push_back(m); 

      pads[ipol]->cd(1)->cd(2); 
      TPaveText * pt  = new TPaveText(i/double(last.nPeaks[ipol]),0,(i+1)/double(last.nPeaks[ipol]),1); 
      delete_list.push_back(pt); 
      pt->AddText(TString::Format("#phi: %0.2f (rough) , %0.3f (fine)", rough_peaks[ipol][i].first, last.peak[ipol][i].phi)); 
      pt->AddText(TString::Format("#theta: %0.2f (rough) , %0.3f (fine)", -rough_peaks[ipol][i].second, last.peak[ipol][i].theta)); 
      pt->AddText(TString::Format("peak val: %f", last.peak[ipol][i].value)); 
      pt->AddText(TString::Format("peak_{hilbert} (coherent): %0.3f", last.coherent[ipol][i].peakHilbert)); 
      pt->AddText(TString::Format("Stokes: (coherent): (%0.3g, %0.3g, %0.3g, %0.3g)", last.coherent[ipol][i].I, last.coherent[ipol][i].Q, last.coherent[ipol][i].U, last.coherent[ipol][i].V));
//      pt->AddText(TString::Format("BW: (coherent): %0.3f", last.coherent[ipol][i].bandwidth)); 
      pt->AddText(TString::Format("position: %0.3f N, %0.3f E, %0.3f m", last.peak[ipol][i].latitude, last.peak[ipol][i].longitude, last.peak[ipol][i].altitude)); 
      pt->Draw(); 
    }


    pads[ipol]->cd(2)->Divide(last.nPeaks[ipol], 3); 

    for (int i = 0; i < last.nPeaks[ipol]; i++) 
    {
      pads[ipol]->cd(2)->cd(i+1); 
      zoomed_correlation_maps[ipol][i]->SetTitle(TString::Format("Zoomed peak %d", i+1)); 
      zoomed_correlation_maps[ipol][i]->Draw("colz"); 
      const AnitaEventSummary::PointingHypothesis & p = last.peak[ipol][i]; 
      TMarker * m = new TMarker(p.phi, -p.theta,2); 
      delete_list.push_back(m); 

      double angle = 90. / TMath::Pi()* atan2(2*p.rho * p.sigma_theta * p.sigma_phi, p.sigma_phi * p.sigma_phi - p.sigma_theta * p.sigma_theta);
      TEllipse *el = new TEllipse(p.phi, -p.theta, p.sigma_phi, p.sigma_theta, 0, 360, angle); 
      delete_list.push_back(el); 

      el->SetFillStyle(0); 
      el->SetLineColor(3); 
      el->Draw(); 
      m->SetMarkerSize(2); 
      m->SetMarkerColor(3); 
      m->Draw(); 

      pads[ipol]->cd(2)->cd(i+last.nPeaks[ipol]+1); 


      ((TGraph*) coherent[ipol][i]->even())->SetTitle(TString::Format ( "Coherent%s (+ xpol) %d", interactive_deconvolved ? " / Deconvolved (+ xpol)" : "", i+1)); 
      coherent[ipol][i]->drawEven("al"); 

      if (interactive_deconvolved)
      {
        deconvolved[ipol][i]->drawEven("lsame"); 
      }

      coherent_xpol[ipol][i]->drawEven("lsame"); 

      if (interactive_xpol_deconvolved)
      {
        deconvolved_xpol[ipol][i]->drawEven("lsame"); 
      }


      pads[ipol]->cd(2)->cd(i+2*last.nPeaks[ipol]+1); 


      (((TGraph*)coherent_power[ipol][i]))->SetTitle(TString::Format ( "Power Coherent%s (+ xpol) %d", interactive_deconvolved ? " / Deconvolved (+ xpol)" : "", i+1)); 
      ((TGraph*)coherent_power[ipol][i])->Draw("al"); 

      if (interactive_deconvolved)
      {
        ((TGraph*)deconvolved_power[ipol][i])->Draw("lsame");; 
      }

      ((TGraph*)coherent_power_xpol[ipol][i])->Draw("lsame"); 

      if (interactive_xpol_deconvolved)
      {
        ((TGraph*)deconvolved_power_xpol[ipol][i])->Draw("lsame"); 
      }
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




