#include "Analyzer.h" 
#include "AnalysisConfig.h" 
#include "TPaveText.h"
#include "AnalysisWaveform.h"
#include "AnitaVersion.h" 
#include "TCanvas.h"
#include "ImpulsivityMeasure.h" 
#include "TStyle.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "FilteredAnitaEvent.h"
#include "UCUtil.h"
#include "DigitalFilter.h" 
#include "FFTtools.h" 
#include "UsefulAdu5Pat.h"
#include "RawAnitaHeader.h" 
#include "PeakFinder.h"
#include "UCFlags.h"
#include "ShapeParameters.h" 
#include "SpectrumParameters.h" 
#include "TF1.h" 
#include "TGraphErrors.h"
#include "simpleStructs.h"
#include "UCImageTools.h"

#ifdef UCORRELATOR_OPENMP
#include <omp.h>
#include "TThread.h" 

#define SECTIONS _Pragma("omp parallel sections")
#define SECTION _Pragma("omp section") 

#else 

#define SECTIONS if(true) 
#define SECTION if(true) 

#endif 

static UCorrelator::AnalysisConfig defaultConfig; 
static int instance_counter = 0; 

#ifndef DEG2RAD
#define DEG2RAD (TMath::Pi()/ 180) 
#endif

#ifndef RAD2DEG
#define RAD2DEG (180 / TMath::Pi()) 
#endif






UCorrelator::Analyzer::Analyzer(const AnalysisConfig * conf, bool interactive_mode) 
  : cfg(conf ? conf: &defaultConfig),
    corr(cfg->correlator_nphi,0,360,  cfg->correlator_ntheta, -cfg->correlator_theta_lowest, cfg->correlator_theta_highest) , 
    responses(AnalysisConfig::getResponseString(cfg->response_option), cfg->response_npad, cfg->deconvolution_method), 
    wfcomb(cfg->combine_nantennas, cfg->combine_npad, true, cfg->response_option!=AnalysisConfig::ResponseNone, &responses), 
    wfcomb_xpol(cfg->combine_nantennas, cfg->combine_npad, true, cfg->response_option!=AnalysisConfig::ResponseNone, &responses), 
    wfcomb_filtered(cfg->combine_nantennas, cfg->combine_npad, false, cfg->response_option!=AnalysisConfig::ResponseNone, &responses), 
    wfcomb_xpol_filtered(cfg->combine_nantennas, cfg->combine_npad, false, cfg->response_option!=AnalysisConfig::ResponseNone, &responses), 
    interactive(interactive_mode)
{
#ifdef UCORRELATOR_OPENMP
  TThread::Initialize(); 
#endif

  zoomed = new TH2D(TString::Format("zoomed_%d", instance_counter), "Zoomed!", cfg->zoomed_nphi, 0 ,1, cfg->zoomed_ntheta, 0, 1);
  zoomed->SetDirectory(0); 

  avg_spectra[0] = 0; 
  avg_spectra[1] = 0; 

  corr.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb_xpol.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb_filtered.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb_xpol_filtered.setGroupDelayFlag(cfg->enable_group_delay); 
  instance_counter++; 
  power_filter = new FFTtools::GaussianFilter(2,3) ; //TODO make this configurable

  if (interactive) 
  {
    correlation_maps[0] = 0;
    correlation_maps[1] = 0;

    for (int i = 0; i < cfg->nmaxima; i++)
    {
      zoomed_correlation_maps[0].push_back(0); 
      zoomed_correlation_maps[1].push_back(0); 

      for (int ipol =0; ipol <2;ipol++)
      {
        for (int ifilt = 0; ifilt <2; ifilt++) 
        {
          coherent[ipol][ifilt].push_back(new AnalysisWaveform); 
          deconvolved[ipol][ifilt].push_back(new AnalysisWaveform); 
          coherent_xpol[ipol][ifilt].push_back(new AnalysisWaveform); 
          deconvolved_xpol[ipol][ifilt].push_back(new AnalysisWaveform); 
          coherent_power[ipol][ifilt].push_back(new TGraphAligned); 
          deconvolved_power[ipol][ifilt].push_back(new TGraphAligned); 
          coherent_power_xpol[ipol][ifilt].push_back(new TGraphAligned); 
          deconvolved_power_xpol[ipol][ifilt].push_back(new TGraphAligned); 
        }
      }
    }
  }
}





void UCorrelator::Analyzer::analyze(const FilteredAnitaEvent * event, AnitaEventSummary * summary, const TruthAnitaEvent * truth) 
{
  
  const RawAnitaHeader * hdr = event->getHeader(); 

  //we need a UsefulAdu5Pat for this event
  UsefulAdu5Pat * pat =  (UsefulAdu5Pat*) event->getGPS();  //unconstifying it .. hopefully that won't cause problems
 
  /* Initialize the summary */ 
  summary = new (summary) AnitaEventSummary(hdr, (UsefulAdu5Pat*) event->getGPS(),truth); 

  //check for saturation
  uint64_t saturated[2] = {0,0}; 
  event->checkSaturation( &saturated[AnitaPol::kHorizontal], 
                          &saturated[AnitaPol::kVertical], 
                          cfg->saturation_threshold); 

  //also disable missing sectors 
  UCorrelator::flags::checkEmpty(event->getUsefulAnitaEvent(), &saturated[AnitaPol::kHorizontal], &saturated[AnitaPol::kVertical]); 

  // loop over wanted polarizations 
  for (int pol = cfg->start_pol; pol <= cfg->end_pol; pol++) 
  {

    UShort_t triggeredPhi = AnitaPol::AnitaPol_t(pol) == AnitaPol::kHorizontal ? event->getHeader()->l3TrigPatternH : event->getHeader()->l3TrigPattern; 
    UShort_t triggeredPhiXpol = AnitaPol::AnitaPol_t(pol) == AnitaPol::kVertical ? event->getHeader()->l3TrigPatternH : event->getHeader()->l3TrigPattern; 


    UShort_t maskedL2 = 0;   
    UShort_t maskedPhi = 0 ; 
    UShort_t maskedL2Xpol = 0;   
    UShort_t maskedPhiXpol = 0 ; 

    if (cfg->use_offline_mask) 
    {
      maskedPhi = event->getHeader()->getPhiMaskOffline(AnitaPol::AnitaPol_t(pol));
      maskedPhiXpol = event->getHeader()->getPhiMaskOffline(AnitaPol::AnitaPol_t(1-pol));

      if (AnitaVersion::get() == 3) 
      {
        maskedL2 = event->getHeader()->getL1MaskOffline(AnitaPol::AnitaPol_t(pol));
        maskedL2Xpol = event->getHeader()->getL1MaskOffline(AnitaPol::AnitaPol_t(1-pol));
      }
      else 
      {

        //for now do this  until we merge anita_3and4 to master
#ifdef MULTIVERSION_ANITA_ENABLED
        maskedL2 = event->getHeader()->getL2Mask(); 
#else
#if VER_EVENT_HEADER>=40
        maskedL2 = event->getHeader()->l2TrigMask;
#else
        maskedL2 = event->getHeader()->l1TrigMask;
#endif
#endif 
        maskedL2Xpol = maskedL2; 

      }
    }
    else
    {
      maskedPhi = AnitaPol::AnitaPol_t(pol) == AnitaPol::kHorizontal ? event->getHeader()->phiTrigMaskH : event->getHeader()->phiTrigMask; 
      maskedPhiXpol = AnitaPol::AnitaPol_t(pol) == AnitaPol::kVertical ? event->getHeader()->phiTrigMaskH : event->getHeader()->phiTrigMask; 

      if (AnitaVersion::get() == 3) 
      {
        maskedL2 = event->getHeader()->getL1Mask(AnitaPol::AnitaPol_t(pol));
        maskedL2Xpol = event->getHeader()->getL1Mask(AnitaPol::AnitaPol_t(1-pol));
      }
      else
      {

        //for now do this  until we merge anita_3and4 to master
#ifdef MULTIVERSION_ANITA_ENABLED
        maskedL2 = event->getHeader()->getL2Mask(); 
#else
#if VER_EVENT_HEADER>=40
        maskedL2 = event->getHeader()->l2TrigMask;
#else
        maskedL2 = event->getHeader()->l1TrigMask;
#endif
#endif
 
        maskedL2Xpol = maskedL2; 
      }


    }


    //alright, we need a little song and dance here to combine Phi masks and L2 masks 
    //
    //  Someone should check my logic here. 
    //
    // An L2 mask effectively masks both phi sector N and N+1 (since it'll prevent L3 triggers in both). 
    // An L3 mask has contributions from both phi sector N and N-1. 
    //
    // I'm going to take the aggressive approach where an L3 mask means we mark a pointing hypothesis as masked
    // if it falls in phi sector N or N-1. 
   

#ifdef __cpp_static_assert
    static_assert(sizeof(maskedPhi == NUM_PHI),"masked phi must be same size as num phi "); 
#endif
     
    maskedPhi |= ( (maskedPhi >> 1) | ( maskedPhi << (NUM_PHI-1))) ; 
    //or with l2 mask
    maskedPhi |= maskedL2; 


#ifdef __cpp_static_assert
    static_assert(sizeof(maskedPhiXpol == NUM_PHI),"masked phi xpol must be same size as num phi "); 
#endif

    //ditto for xpol 
    maskedPhiXpol |=  (maskedPhiXpol >> 1) | ((maskedPhiXpol << (NUM_PHI-1))); 
    maskedPhiXpol |= maskedL2Xpol; 

    //TODO: check this 
    TVector2 triggerAngle(0,0); 

    int ntriggered = __builtin_popcount(triggeredPhi); 
    int ntriggered_xpol = __builtin_popcount(triggeredPhi); 

    UShort_t which_trigger = ntriggered ? triggeredPhi : triggeredPhiXpol; 
    int which_ntriggered = ntriggered ?: ntriggered_xpol; 

    const double phi_sector_width = 360. / NUM_PHI; 

    for (int i = 0; i < NUM_PHI; i++) 
    {
      if (which_trigger & (1 << i))
      {
        //TODO: this 45 is hardcoded here. Should come from GeomTool or something... 
         double ang = (i * phi_sector_width - 45) * TMath::Pi()/180;
         triggerAngle += TVector2(cos(ang), sin(ang)) / which_ntriggered; 
      }
    }

    double avgHwAngle = triggerAngle.Phi() * RAD2DEG; 

    // tell the correlator not to use saturated events and make the correlation map
    corr.setDisallowedAntennas(saturated[pol]); 
    corr.compute(event, AnitaPol::AnitaPol_t(pol)); 

    //compute RMS of correlation map 
    //    maprms = corr.getHist()->GetRMS(3); //This doesn't work!  Probably because ROOT is dumb
    maprms = UCorrelator::getZRMS(corr.getHist());


    // Find the isolated peaks in the image 
    peakfinder::RoughMaximum maxima[cfg->nmaxima]; 
    int npeaks = UCorrelator::peakfinder::findIsolatedMaxima((const TH2D*) corr.getHist(), cfg->peak_isolation_requirement, cfg->nmaxima, maxima, cfg->use_bin_center); 
//    printf("npeaks: %d\n", npeaks); 
    summary->nPeaks[pol] = npeaks; 

    rough_peaks[pol].clear(); 


    //get the average spectra 
    if (!avg_spectra[pol]) 
    {
      avg_spectra[pol] = new TGraph; 
      avg_spectra[pol]->GetXaxis()->SetTitle("Frequency (GHz)"); 
      avg_spectra[pol]->GetYaxis()->SetTitle("Power (dBish)"); 
      avg_spectra[pol]->SetTitle(TString::Format("Average spectra for %s", pol == AnitaPol::kHorizontal ? "HPol" : "VPol")); 
    }

    event->getMedianSpectrum(avg_spectra[pol], AnitaPol::AnitaPol_t(pol),0.5); 

    // Loop over found peaks 
    for (int i = 0; i < npeaks; i++) 
    {
      // zoom in on the values 
//      printf("rough phi:%f, rough theta: %f\n", maxima[i].x, -maxima[i].y); 


      fillPointingInfo(maxima[i].x, maxima[i].y, &summary->peak[pol][i], pat, avgHwAngle, triggeredPhi, maskedPhi, triggeredPhiXpol, maskedPhiXpol); 

      if (interactive) 
      {
        if (zoomed_correlation_maps[pol][i]) delete zoomed_correlation_maps[pol][i]; 
        zoomed_correlation_maps[pol][i] = new gui::Map(*zoomed, event, &wfcomb, &wfcomb_filtered,AnitaPol::AnitaPol_t(pol), summary); 
        zoomed_correlation_maps[pol][i]->SetName(TString::Format("zoomed_%d_%d", pol,i)); 
      }


      //fill in separation 
      summary->peak[pol][i].phi_separation = 1000; 
      for (int j = 0; j < i; j++)
      {
        summary->peak[pol][i].phi_separation = TMath::Min(summary->peak[pol][i].phi_separation, fabs(FFTtools::wrap(summary->peak[pol][i].phi - summary->peak[pol][j].phi, 360, 0))); 
      }

//      printf("phi:%f, theta:%f\n", summary->peak[pol][i].phi, summary->peak[pol][i].theta); 


    }

    for (int i = 0; i < npeaks; i++) 
    {
      rough_peaks[pol].push_back(std::pair<double,double>(maxima[i].x, maxima[i].y)); 
      //now make the combined waveforms 
     
      
SECTIONS
{
SECTION
      wfcomb.combine(summary->peak[pol][i].phi, summary->peak[pol][i].theta, event, (AnitaPol::AnitaPol_t) pol, saturated[pol]); 
SECTION
      wfcomb_xpol.combine(summary->peak[pol][i].phi, summary->peak[pol][i].theta, event, (AnitaPol::AnitaPol_t) (1-pol), saturated[pol]); 
SECTION
      wfcomb_filtered.combine(summary->peak[pol][i].phi, summary->peak[pol][i].theta, event, (AnitaPol::AnitaPol_t) pol, saturated[pol]); 
SECTION
      wfcomb_xpol_filtered.combine(summary->peak[pol][i].phi, summary->peak[pol][i].theta, event, (AnitaPol::AnitaPol_t) (1-pol), saturated[pol]); 
}

SECTIONS 
{
SECTION
      fillWaveformInfo(wfcomb.getCoherent(), wfcomb_xpol.getCoherent(), wfcomb.getCoherentAvgSpectrum(), &summary->coherent[pol][i], (AnitaPol::AnitaPol_t) pol); 
SECTION
      fillWaveformInfo(wfcomb.getDeconvolved(), wfcomb_xpol.getDeconvolved(), wfcomb.getDeconvolvedAvgSpectrum(), &summary->deconvolved[pol][i],  (AnitaPol::AnitaPol_t)pol); 
SECTION
      fillWaveformInfo(wfcomb_filtered.getCoherent(), wfcomb_xpol_filtered.getCoherent(), wfcomb_filtered.getCoherentAvgSpectrum(), &summary->coherent_filtered[pol][i], (AnitaPol::AnitaPol_t) pol); 
SECTION
      fillWaveformInfo(wfcomb_filtered.getDeconvolved(), wfcomb_xpol_filtered.getDeconvolved(), wfcomb_filtered.getDeconvolvedAvgSpectrum(), &summary->deconvolved_filtered[pol][i],  (AnitaPol::AnitaPol_t)pol); 
}



      if (interactive) //copy everything
      {
        coherent[pol][0][i]->~AnalysisWaveform(); 
        coherent[pol][0][i] = new (coherent[pol][0][i]) AnalysisWaveform(*wfcomb.getCoherent()); 

        coherent[pol][1][i]->~AnalysisWaveform(); 
        coherent[pol][1][i] = new (coherent[pol][1][i]) AnalysisWaveform(*wfcomb_filtered.getCoherent()); 

        coherent_power[pol][0][i]->~TGraphAligned(); 
        coherent_power[pol][0][i] = new (coherent_power[pol][0][i]) TGraphAligned(*wfcomb.getCoherentAvgSpectrum()); 
        coherent_power[pol][0][i]->dBize(); 

        coherent_power[pol][1][i]->~TGraphAligned(); 
        coherent_power[pol][1][i] = new (coherent_power[pol][1][i]) TGraphAligned(*wfcomb_filtered.getCoherentAvgSpectrum()); 
        coherent_power[pol][1][i]->dBize(); 


        if (wfcomb.getDeconvolved())
        {
          deconvolved[pol][0][i]->~AnalysisWaveform(); 
          deconvolved[pol][0][i] = new (deconvolved[pol][0][i]) AnalysisWaveform(*wfcomb.getDeconvolved()); 
          deconvolved[pol][0][i]->updateEven()->SetLineColor(2); 
          deconvolved[pol][0][i]->updateEven()->SetMarkerColor(2); 

          deconvolved[pol][1][i]->~AnalysisWaveform(); 
          deconvolved[pol][1][i] = new (deconvolved[pol][1][i]) AnalysisWaveform(*wfcomb_filtered.getDeconvolved()); 
          deconvolved[pol][1][i]->updateEven()->SetLineColor(2); 
          deconvolved[pol][1][i]->updateEven()->SetMarkerColor(2); 

          deconvolved_power[pol][0][i]->~TGraphAligned(); 
          deconvolved_power[pol][0][i] = new (deconvolved_power[pol][0][i]) TGraphAligned(*wfcomb.getDeconvolvedAvgSpectrum()); 
          deconvolved_power[pol][0][i]->dBize(); 
          deconvolved_power[pol][0][i]->SetLineColor(2); 

          deconvolved_power[pol][1][i]->~TGraphAligned(); 
          deconvolved_power[pol][1][i] = new (deconvolved_power[pol][1][i]) TGraphAligned(*wfcomb_filtered.getDeconvolvedAvgSpectrum()); 
          deconvolved_power[pol][1][i]->dBize(); 
          deconvolved_power[pol][1][i]->SetLineColor(2); 
 
          interactive_deconvolved = true; 
        }
        else
        {
          interactive_deconvolved = false; 
        }


        coherent_xpol[pol][0][i]->~AnalysisWaveform(); 
        coherent_xpol[pol][0][i] = new (coherent_xpol[pol][0][i]) AnalysisWaveform(*wfcomb_xpol.getCoherent()); 
        coherent_xpol[pol][0][i]->updateEven()->SetLineColor(11); 
        coherent_xpol[pol][0][i]->updateEven()->SetLineStyle(3); 

        coherent_xpol[pol][1][i]->~AnalysisWaveform(); 
        coherent_xpol[pol][1][i] = new (coherent_xpol[pol][1][i]) AnalysisWaveform(*wfcomb_xpol_filtered.getCoherent()); 
        coherent_xpol[pol][1][i]->updateEven()->SetLineColor(11); 
        coherent_xpol[pol][1][i]->updateEven()->SetLineStyle(3); 
 
        
        coherent_power_xpol[pol][0][i]->~TGraphAligned(); 
        coherent_power_xpol[pol][0][i] = new (coherent_power_xpol[pol][0][i]) TGraphAligned(*wfcomb_xpol.getCoherentAvgSpectrum()); 
        coherent_power_xpol[pol][0][i]->dBize(); 
        coherent_power_xpol[pol][0][i]->SetLineStyle(3); 
        coherent_power_xpol[pol][0][i]->SetLineColor(11); 

        coherent_power_xpol[pol][1][i]->~TGraphAligned(); 
        coherent_power_xpol[pol][1][i] = new (coherent_power_xpol[pol][1][i]) TGraphAligned(*wfcomb_xpol_filtered.getCoherentAvgSpectrum()); 
        coherent_power_xpol[pol][1][i]->dBize(); 
        coherent_power_xpol[pol][1][i]->SetLineStyle(3); 
        coherent_power_xpol[pol][1][i]->SetLineColor(11); 


        if (wfcomb_xpol.getDeconvolved())
        {
          deconvolved_xpol[pol][0][i]->~AnalysisWaveform(); 
          deconvolved_xpol[pol][0][i] = new (deconvolved_xpol[pol][0][i]) AnalysisWaveform(*wfcomb_xpol.getDeconvolved()); 
          deconvolved_xpol[pol][0][i]->updateEven()->SetLineColor(45); 
          deconvolved_xpol[pol][0][i]->updateEven()->SetMarkerColor(45); 
          deconvolved_xpol[pol][0][i]->updateEven()->SetLineStyle(3); 

          deconvolved_power_xpol[pol][0][i]->~TGraphAligned(); 
          deconvolved_power_xpol[pol][0][i] = new (deconvolved_power_xpol[pol][0][i]) TGraphAligned(*wfcomb_xpol.getDeconvolvedAvgSpectrum()); 
          deconvolved_power_xpol[pol][0][i]->dBize(); 
          deconvolved_power_xpol[pol][0][i]->SetLineStyle(3); 
          deconvolved_power_xpol[pol][0][i]->SetLineColor(46); 

          deconvolved_xpol[pol][1][i]->~AnalysisWaveform(); 
          deconvolved_xpol[pol][1][i] = new (deconvolved_xpol[pol][1][i]) AnalysisWaveform(*wfcomb_xpol_filtered.getDeconvolved()); 
          deconvolved_xpol[pol][1][i]->updateEven()->SetLineColor(45); 
          deconvolved_xpol[pol][1][i]->updateEven()->SetMarkerColor(45); 
          deconvolved_xpol[pol][1][i]->updateEven()->SetLineStyle(3); 

          deconvolved_power_xpol[pol][1][i]->~TGraphAligned(); 
          deconvolved_power_xpol[pol][1][i] = new (deconvolved_power_xpol[pol][1][i]) TGraphAligned(*wfcomb_xpol_filtered.getDeconvolvedAvgSpectrum()); 
          deconvolved_power_xpol[pol][1][i]->dBize(); 
          deconvolved_power_xpol[pol][1][i]->SetLineStyle(3); 
          deconvolved_power_xpol[pol][1][i]->SetLineColor(46); 



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
       if (correlation_maps[pol]) delete correlation_maps[pol];
       correlation_maps[pol] = new gui::Map(*corr.getHist(), event, &wfcomb, &wfcomb_filtered,AnitaPol::AnitaPol_t(pol), summary ); 
    }
  }

  fillFlags(event, &summary->flags, pat); 

  if (truth)
  { 
//    SECTIONS
    {
//      SECTION
      wfcomb.combine(summary->mc.phi, summary->mc.theta, event, AnitaPol::kHorizontal, 0); 
//      SECTION 
      wfcomb_xpol.combine(summary->mc.phi, summary->mc.theta, event, AnitaPol::kVertical, 0); 
    }

//    SECTIONS
    {
//      SECTION
      fillWaveformInfo(wfcomb.getCoherent(), wfcomb_xpol.getCoherent(), wfcomb.getCoherentAvgSpectrum(), &(summary->mc.wf[AnitaPol::kHorizontal]), AnitaPol::kHorizontal); 
//      SECTION
      fillWaveformInfo(wfcomb_xpol.getCoherent(), wfcomb.getCoherent(), wfcomb_xpol.getCoherentAvgSpectrum(), &(summary->mc.wf[AnitaPol::kVertical]), AnitaPol::kVertical); 
    }
  }

  if (interactive) last = *summary; 

}

static bool outside(const TH2 * h, double x, double y) 
{

  return x > h->GetXaxis()->GetXmax() || 
         x < h->GetXaxis()->GetXmin() || 
         y < h->GetYaxis()->GetXmin() ||
         y > h->GetYaxis()->GetXmax(); 

}

void UCorrelator::Analyzer::fillPointingInfo(double rough_phi, double rough_theta, AnitaEventSummary::PointingHypothesis * point, 
                                             UsefulAdu5Pat * pat, double hwPeakAngle, UShort_t triggered_sectors, UShort_t masked_sectors, UShort_t triggered_sectors_xpol, UShort_t masked_sectors_xpol)
{
      corr.computeZoomed(rough_phi, rough_theta, cfg->zoomed_nphi, cfg->zoomed_dphi,  cfg->zoomed_ntheta, cfg->zoomed_dtheta, cfg->zoomed_nant, zoomed); 

      //get pointer to the pointing hypothesis we are about to fill 

      // This will fill in phi, theta, value, var_theta, var_phi and covar 
      
      peakfinder::FineMaximum max; 
      switch (cfg->fine_peak_finding_option)
      {
        case AnalysisConfig::FinePeakFindingAbby: 
          UCorrelator::peakfinder::doInterpolationPeakFindingAbby(zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingBicubic: 
          UCorrelator::peakfinder::doInterpolationPeakFindingBicubic(zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingHistogram: 
          UCorrelator::peakfinder::doPeakFindingHistogram(zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingQuadraticFit16: 
          UCorrelator::peakfinder::doPeakFindingQuadratic16(zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingQuadraticFit25: 
          UCorrelator::peakfinder::doPeakFindingQuadratic25(zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingGaussianFit: 
          UCorrelator::peakfinder::doPeakFindingGaussian(zoomed, &max); 
          break; 
        case AnalysisConfig::FinePeakFindingQuadraticFit9: 
        default: 
          UCorrelator::peakfinder::doPeakFindingQuadratic9(zoomed, &max); 
          break; 
      }; 


      //Check to make sure that fine max isn't OUTSIDE of zoomed window
      // If it is, revert to very stupid method  of just using histogram 
      
      if (outside(zoomed, max.x, max.y))
      {
        UCorrelator::peakfinder::doPeakFindingHistogram(zoomed, &max); 
      }
      


      
      max.copyToPointingHypothesis(point); 

      //snr is ratio of point value to map rms
      point->snr = point->value / maprms; 
      point->dphi_rough = FFTtools::wrap(point->phi - rough_phi, 360,0); 
      point->dtheta_rough = FFTtools::wrap(point->theta - (-rough_theta), 360,0); //sign reversal. doh. 

      point->hwAngle = FFTtools::wrap(point->phi - hwPeakAngle,360,0); 

      //TODO: I don't believe this really yet
      int sector = 2+fmod(point->phi + 11.25,360) / 22.5; 

      point->masked = masked_sectors & ( 1 << sector); 
      point->triggered = triggered_sectors & ( 1 << sector); 
      point->masked_xpol = masked_sectors_xpol & ( 1 << sector); 
      point->triggered_xpol = triggered_sectors_xpol & ( 1 << sector); 

      //Compute intersection with continent, or set values to -9999 if no intersection
      if (!pat->traceBackToContinent(point->phi * DEG2RAD, point->theta * DEG2RAD, &point->longitude, &point->latitude, &point->altitude, &point->theta_adjustment_needed)) 
      {
        point->latitude = -9999; 
        point->longitude = -9999;  
        point->altitude = -9999; 
        point->distanceToSource = -9999; 
        point->theta_adjustment_needed = -9999; 
      }
      else
      {
        point->distanceToSource=pat->getDistanceFromSource(point->latitude, point->longitude, point->altitude); 
        point->theta_adjustment_needed *= RAD2DEG; 
      }
}


void UCorrelator::Analyzer::fillWaveformInfo(const AnalysisWaveform * wf, const AnalysisWaveform * xpol_wf, const TGraph* pwr, AnitaEventSummary::WaveformInfo * info, AnitaPol::AnitaPol_t pol)
{
  if (!wf || wf->Neven() == 0)
  {
    if (wf && !wf->Neven()) 
      fprintf(stderr,"wf passed to fillWaveformInfo has no points\n");  

    memset(info, 0, sizeof(AnitaEventSummary::WaveformInfo)); 
    return; 
  }
  const TGraphAligned * even = wf->even(); 
  const TGraphAligned * xpol_even= xpol_wf->even(); 
  int peakBin;

  int peakHilbertBin; 
  info->peakVal = FFTtools::getPeakVal((TGraph*) even,&peakBin); 
  info->xPolPeakVal = FFTtools::getPeakVal( xpol_even); 
  info->peakHilbert = FFTtools::getPeakVal((TGraph*) wf->hilbertEnvelope(),&peakHilbertBin); 
  double minHilbert = *std::min_element(wf->hilbertEnvelope()->GetY(), wf->hilbertEnvelope()->GetY() + wf->Neven()); 

  info->xPolPeakHilbert = FFTtools::getPeakVal((TGraph*) xpol_wf->hilbertEnvelope()); 
  info->numAntennasInCoherent = cfg->combine_nantennas; 

  info->totalPower = even->getSumV2(); 
  info->totalPowerXpol = xpol_even->getSumV2(); 
  info->peakTime = even->GetX()[peakHilbertBin]; 

  double hilbertRange = info->peakHilbert - minHilbert; 

  info->riseTime_10_90 = shape::getRiseTime((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, minHilbert + 0.9*hilbertRange,peakHilbertBin); 
  info->riseTime_10_50 = shape::getRiseTime((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, minHilbert + 0.5*hilbertRange,peakHilbertBin); 
  info->fallTime_90_10 = shape::getFallTime((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, minHilbert + 0.9*hilbertRange,peakHilbertBin); 
  info->fallTime_50_10 = shape::getFallTime((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, minHilbert + 0.5*hilbertRange,peakHilbertBin); 

  int ifirst, ilast; 
  info->width_50_50 = shape::getWidth((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.5*hilbertRange, &ifirst, &ilast,peakHilbertBin); 
  info->power_50_50 = even->getSumV2(ifirst, ilast); 
  even->getMoments(sizeof(info->peakMoments)/sizeof(double), info->peakTime, info->peakMoments); 
  info->width_10_10 = shape::getWidth((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, &ifirst, &ilast,peakHilbertBin); 
  info->power_10_10 = even->getSumV2(ifirst, ilast); 


  if (ifirst < 0) ifirst = 0; 
  if (ilast < 0) ilast = wf->Neven()-1; 
  int nstokes = ilast-ifirst+1 ; 

  if (pol == AnitaPol::kHorizontal)
  {

    FFTtools::stokesParameters(nstokes,
                               even->GetY()+ifirst, 
                               wf->hilbertTransform()->even()->GetY()+ifirst, 
                               xpol_even->GetY()+ifirst, 
                               xpol_wf->hilbertTransform()->even()->GetY()+ifirst, 
                               &(info->I), &(info->Q), &(info->U), &(info->V)); 
  }
  else
  {
    FFTtools::stokesParameters(nstokes, 
                               xpol_even->GetY()+ifirst, 
                               xpol_wf->hilbertTransform()->even()->GetY()+ifirst, 
                               even->GetY()+ifirst, 
                               wf->hilbertTransform()->even()->GetY()+ifirst, 
                               &(info->I), &(info->Q), &(info->U), &(info->V)); 
 
  }

  info->impulsivityMeasure = impulsivity::impulsivityMeasure(wf); 

  double dt = wf->deltaT(); 
  double t0 = even->GetX()[0]; 

  int i0 = TMath::Max(0.,floor((cfg->noise_estimate_t0 - t0)/dt)); 
  int i1 = TMath::Min(even->GetN()-1.,ceil((cfg->noise_estimate_t1 - t0)/dt)); 
  int n = i1 - i0 + 1; 
//  printf("%d-%d -> %d \n", i0, i1, n); 

  double rms = TMath::RMS(n, even->GetY() + i0); 
  
  info->snr = info->peakVal / rms;
  TGraphAligned power(pwr->GetN(),pwr->GetX(),pwr->GetY()); 
  power.dBize(); 

  if (power_filter)
  {
    power_filter->filterGraph(&power); 
  }

  spectrum::fillSpectrumParameters(&power, avg_spectra[pol], info, cfg); 
}


UCorrelator::Analyzer::~Analyzer()
{

  delete zoomed; 
  if (interactive)
  {
    delete correlation_maps[0];
    delete correlation_maps[1]; 


    for (int pol = 0; pol < 2; pol++)
    {
      for (int i = 0; i < cfg->nmaxima; i++)
      {
          delete zoomed_correlation_maps[pol][i];

        for (int ifilt = 0; ifilt < 2; ifilt++) 
        {

          delete coherent[pol][ifilt][i];
          delete deconvolved[pol][ifilt][i];

          delete coherent_xpol[pol][ifilt][i];
          delete deconvolved_xpol[pol][ifilt][i];

          delete coherent_power[pol][ifilt][i];
          delete deconvolved_power[pol][ifilt][i];

          delete coherent_power_xpol[pol][ifilt][i];
          delete deconvolved_power_xpol[pol][ifilt][i];
        }

      }
    }

    clearInteractiveMemory(1); 
  }

  if (power_filter)
    delete power_filter; 
}

void UCorrelator::Analyzer::clearInteractiveMemory(double frac) const
{

  for (unsigned i = (1-frac) * delete_list.size(); i < delete_list.size(); i++) 
  {
    delete delete_list[i]; 
  }

  delete_list.clear(); 
}


/* Nevermind this... wanted to zoom in on analyzer canvas on click, but too much work :( 
static void setOnClickHandler(TPad * pad) 
{

}
*/





void UCorrelator::Analyzer::drawSummary(TPad * ch, TPad * cv) const
{
  TPad * pads[2] = {ch,cv}; 

  clearInteractiveMemory(); 

//  gStyle->SetOptStat(0); 
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
    correlation_maps[ipol]->addRough(rough_peaks[ipol]); 
    correlation_maps[ipol]->Draw("colz"); 

    for (int i = 0; i < last.nPeaks[ipol]; i++) 
    {
      pads[ipol]->cd(1)->cd(2); 
      UCorrelator::gui::SummaryText * pt  = new gui::SummaryText(i, AnitaPol::AnitaPol_t(ipol), this); 
      delete_list.push_back(pt); 
      pt->Draw(); 
    }


    pads[ipol]->cd(2)->Divide(last.nPeaks[ipol], interactive_deconvolved ? 5 : 3); 

    for (int i = 0; i < last.nPeaks[ipol]; i++) 
    {
      pads[ipol]->cd(2)->cd(i+1); 

      zoomed_correlation_maps[ipol][i]->SetTitle(TString::Format("Zoomed peak %d", i+1)); 
      zoomed_correlation_maps[ipol][i]->addFine(last.peak[ipol][i]); 
      zoomed_correlation_maps[ipol][i]->Draw("colz"); 

      pads[ipol]->cd(2)->cd(i+last.nPeaks[ipol]+1); 

      ((TGraph*) coherent[ipol][0][i]->even())->SetTitle(TString::Format ( "Coherent (+ xpol) %d", i+1)); 
      coherent[ipol][0][i]->drawEven("al"); 


      coherent_xpol[ipol][0][i]->drawEven("lsame"); 


      pads[ipol]->cd(2)->cd(i+2*last.nPeaks[ipol]+1); 


      (((TGraph*)coherent_power[ipol][0][i]))->SetTitle(TString::Format ( "Power Coherent (+ xpol) %d", i+1)); 
      ((TGraph*)coherent_power[ipol][0][i])->Draw("al"); 


      ((TGraph*)avg_spectra[ipol])->SetLineColor(2); 
      ((TGraph*)avg_spectra[ipol])->Draw("lsame"); 

      /*
      TF1 * spectral_slope = new TF1(TString::Format("__slope_%d", i), "pol1",0.2,0.7); 
      spectral_slope->SetParameter(0, last.coherent[ipol][i].spectrumIntercept) ;
      spectral_slope->SetParameter(1, last.coherent[ipol][i].spectrumSlope) ;


      TGraphErrors *gbw = new TGraphErrors(AnitaEventSummary::peaksPerSpectrum); 
      gbw->SetTitle("Bandwidth Peaks"); 
      for (int bwpeak = 0; bwpeak < AnitaEventSummary::peaksPerSpectrum; bwpeak++) 
      {
        double bwf = last.coherent[ipol][i].peakFrequency[bwpeak]; 
        gbw->SetPoint(bwpeak, bwf, avg_spectra[ipol]->Eval(bwf)+ last.coherent[ipol][i].peakPower[bwpeak]); 
        gbw->SetPointError(bwpeak  , last.coherent[ipol][i].bandwidth[bwpeak]/2,0);
      }
      gbw->SetMarkerColor(4); 
      gbw->SetMarkerStyle(20); 
      gbw->Draw("psame"); 


      delete_list.push_back(spectral_slope); 
      delete_list.push_back(gbw);

      */  
      
      if (interactive_deconvolved)
      {
        pads[ipol]->cd(2)->cd(i+3*last.nPeaks[ipol]+1); 
        ((TGraph*) deconvolved[ipol][0][i]->even())->SetTitle(TString::Format ( "Deconvolved (+ xpol) %d", i+1)); 
        deconvolved[ipol][0][i]->drawEven("alp"); 
        if (interactive_xpol_deconvolved)
        {
            deconvolved_xpol[ipol][0][i]->drawEven("lsame"); 
        }

        pads[ipol]->cd(2)->cd(i+4*last.nPeaks[ipol]+1); 

        (((TGraph*)deconvolved_power[ipol][0][i]))->SetTitle(TString::Format ( "Power Deconvolved (+ xpol) %d", i+1)); 
        ((TGraph*)deconvolved_power[ipol][0][i])->Draw();; 
        if (interactive_xpol_deconvolved)
        {
          ((TGraph*)deconvolved_power_xpol[ipol][0][i])->Draw("lsame"); 
        }


      }
       
    }
  }
}

void UCorrelator::Analyzer::fillFlags(const FilteredAnitaEvent * fae, AnitaEventSummary::EventFlags * flags, UsefulAdu5Pat * pat) 
{

  flags->nadirFlag = true; // we should get rid of htis I guess? 

  
  flags->meanPower[0] = fae->getAveragePower(); 
  flags->medianPower[0] = fae->getMedianPower(); 
  flags->meanPowerFiltered[0] = fae->getAveragePower(AnitaPol::kNotAPol, AnitaRing::kNotARing, true); 
  flags->medianPowerFiltered[0] = fae->getMedianPower(AnitaPol::kNotAPol, AnitaRing::kNotARing, true); 

  for (int ring = 0; ring <AnitaRing::kNotARing; ring++)
  {
    flags->meanPower[1+ring] = fae->getAveragePower(AnitaPol::kNotAPol, AnitaRing::AnitaRing_t(ring)); 
    flags->medianPower[1+ring] = fae->getMedianPower(AnitaPol::kNotAPol, AnitaRing::AnitaRing_t(ring)); 
    flags->meanPowerFiltered[1+ring] = fae->getAveragePower(AnitaPol::kNotAPol, AnitaRing::AnitaRing_t(ring), true); 
    flags->medianPowerFiltered[1+ring] = fae->getMedianPower(AnitaPol::kNotAPol, AnitaRing::AnitaRing_t(ring), true); 
  }


  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {
     fae->getMinMaxRatio(AnitaPol::AnitaPol_t(pol), &flags->maxBottomToTopRatio[pol], &flags->minBottomToTopRatio[pol], &flags->maxBottomToTopRatioSector[pol], &flags->minBottomToTopRatioSector[pol], AnitaRing::kBottomRing, AnitaRing::kTopRing,1); 
  }


  if ( isLDB(fae->getHeader(), cfg))
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

  // more than 80 percent filterd out 
  flags->strongCWFlag = flags->medianPowerFiltered[0] / flags->medianPower[0] < 0.2; 

  flags->isPayloadBlast =  
    (cfg->max_mean_power_filtered && flags->meanPowerFiltered[0] > cfg->max_mean_power_filtered) ||
    (cfg->max_median_power_filtered && flags->medianPowerFiltered[0] > cfg->max_median_power_filtered) ||
    (cfg->max_bottom_to_top_ratio && flags->maxBottomToTopRatio[0] > cfg->max_bottom_to_top_ratio) || 
    (cfg->max_bottom_to_top_ratio && flags->maxBottomToTopRatio[1] > cfg->max_bottom_to_top_ratio); 

  flags->isVarner = false; 
  flags->isVarner2 = false; 

  flags->isGood = !flags->isVarner && !flags->isVarner2 && !flags->strongCWFlag; 

}




