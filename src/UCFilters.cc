#include "FFTtools.h" 
#include "FilterStrategy.h"
#include "UCFilters.h" 
#include "BasicFilters.h" 
#include "TGraph.h" 
#include "AnalysisWaveform.h" 
#include "RawAnitaHeader.h"
#include "TRandom.h"
#include "Adu5Pat.h"
#include "DigitalFilter.h"
#include "AntennaPositions.h" 
#ifdef UCORRELATOR_OPENMP
#include "omp.h"
#endif

void UCorrelator::applyAbbysFilterStrategy(FilterStrategy * strategy) 
{

  //passband between 200 and 1200 MHz 
  strategy->addOperation(new SimplePassBandFilter(0.2, 1.2)); 

  //alfa filter
  strategy->addOperation(new ALFAFilter); 

  // satellite filter if north facing 

  strategy->addOperation(new ConditionalFilterOperation
                             (  new ComplicatedNotchFilter(0.26-0.026, 0.26 + 0.026), 
                                    &UCorrelator::antennaIsNorthFacing, "_north", "if is facing north",true
                             )
                        );  


  //add adaptive filter 
  strategy->addOperation(new AdaptiveFilter(2, getenv("UCORRELATOR_BASELINE_DIR"))); 
}



void UCorrelator::ComplicatedNotchFilter::processOne(AnalysisWaveform* g) 
{
//  printf("ComplicatedNotchFilter::processOne!\n"); 
  int nfreq = g->Nfreq(); 
  double df = g->deltaF();  // this is ordinarily Ghz 

  FFTWComplex * vals = g->updateFreq(); 
  const double kB = 1.38e-23; 
  const double ohms = 50; 

  /* TODO: check this ... factor of a few different from Abby's version */ 
  const double meanpower= 4 * kB * temperature * ohms * (df * 1e9) * TMath::Power(10,gain/10);; 

  double sigma = sqrt(meanpower); 

  for (int i = 0; i < nfreq; i++) 
  {
    if (i * df >= min &&  i * df <= max) 
    {
      double phase = gRandom->Uniform(0, 2*TMath::Pi()); 
      double mag = sigma * sqrt(-2 * log (gRandom->Uniform(0,1))); 
      vals[i].re = mag * cos(phase); 
      vals[i].im = mag * sin(phase); 
    }
  }
}



bool UCorrelator::antennaIsNorthFacing(FilteredAnitaEvent * ev, int ant, AnitaPol::AnitaPol_t pol) 
{
  
  double headingAnt = ev->getGPS()->heading - AntennaPositions::instance()->phiAnt[pol][ant]; 
  if (headingAnt >= 360) headingAnt -=360; 
  if (headingAnt < 0) headingAnt +=360; 
  return (headingAnt  <=120 || headingAnt >= 240); 
}


inline static int getBin(double f, const TGraph * g) 
{
  int bin =   (f- g->GetX()[0]) / (g->GetX()[1] - g->GetX()[0]); 
  if (bin < 0) return 0; 
  if (bin >= g->GetN()) return g->GetN()-1; 
  return bin; 
}


static double getMeanBetween(const TGraph * g, double min, double max) 
{
  double sum = 0; 
  int n = 0; 
  for (int i = getBin(min,g); i <= getBin(max,g); i++) 
  {
      sum += g->GetY()[i]; 
      n++; 
  }

  return sum/n;
}

/*
static void undBize(const TGraph * g) 
{
  for (int i = 0; i < g->GetN(); i++) 
  {
     if (g->GetY()[i] <= -1000) 
     {
       g->GetY()[i] = 0; 
     }
     else
     {
       g->GetY()[i] = TMath::Power(g->GetY()[i]/10,10); 
     }
   }
}
*/

static void dBize(TGraph * g) 
{
  for (int i = 0; i < g->GetN(); i++) 
  {
     if (g->GetY()[i] <=0) 
     {
       g->GetY()[i] = -1000; 
     }
     else
     {
       g->GetY()[i] = 10 * TMath::Log10(g->GetY()[i]); 
     }
   }
}

static void interpolateFreqsBetween(TGraph * g, double fmin, double fmax) 
{

  int bin_before = -1; 
  int bin_after = -1; 

  for (int i = 1; i < g->GetN(); i++) 
  {
    if (bin_before  < 0 && g->GetX()[i] > fmin) 
    {
      bin_before = i-1; 
    }

    if (g->GetX()[i] > fmax) 
    {
      bin_after = i; 
      break;
    }
  }

  double y0 = g->GetY()[bin_before];
  double y1 = g->GetY()[bin_after];

  int diff = (bin_after - bin_before); 


  for (int i = bin_before +1; i < bin_after; i++) 
  {
    g->GetY()[i] = (y0 * i + (diff - i) * y1) / diff; 
  }
}

static std::vector<TString> hpol_freq_names; 
static std::vector<TString> vpol_freq_names; 

const char * UCorrelator::AdaptiveFilter::outputName(unsigned i) const
{

  switch(i) 
  {
    case AnitaPol::kHorizontal: 
      return "mean_freq_hpol";
    case AnitaPol::kVertical: 
      return "mean_freq_vpol";
    case 2+AnitaPol::kHorizontal: 
      return "strongest_cw_hpol"; 
    case 2+AnitaPol::kVertical: 
      return "strongest_cw_vpol"; 
    case 5: 
      return "hpol_freqs"; 
    case 6: 
      return "vpol_freqs"; 
    default: 
      return 0; 
  }

  return 0; 
}


void UCorrelator::AdaptiveFilter::fillOutput(unsigned  i, double *vals) const
{

  switch(i) 
  {
    case 0: 
    case 1: 
      *vals=mean_freq[i]; 
      break; 
    case 2: 
    case 3: 
      *vals=strongest_cw[i-2]; 
      break; 
    case 4: 
      memcpy(vals, freqs[0], nfreq*sizeof(double)); 
      break;
    case 5: 
      memcpy(vals, freqs[1], nfreq*sizeof(double)); 
      break; 
    default: 
      return; 
  }
}

void UCorrelator::AdaptiveFilter::process(FilteredAnitaEvent * event) 
{
//  printf("AdaptiveFilter::process!\n"); 
  /* Here we do all the same processing to the baseline as Abby does */ 
  if (run < 0 || run != event->getHeader()->run)
  {
    run = event->getHeader()->run; 
    if (baseline) delete baseline; 
    baseline = new Baseline(run, navg, baseline_dir); 

    if (hpol_avg) delete hpol_avg; 
    if (vpol_avg) delete vpol_avg; 

    hpol_avg  = new TGraph(*baseline->getBaselineAverageHPol()); 
    vpol_avg = new TGraph(*baseline->getBaselineAverageVPol())  ; 

    //convert to dB 
    dBize(hpol_avg); 
    dBize(vpol_avg); 

    //smooth with box filter
    FFTtools::BoxFilter boxf(7);  
    boxf.filterGraph(hpol_avg); 
    boxf.filterGraph(vpol_avg); 


    //linearly interpolate over satellite hump
    interpolateFreqsBetween(hpol_avg, 0.23, 0.31); 
    interpolateFreqsBetween(vpol_avg, 0.23, 0.31); 


    //means used for normalization
  }



  //find average hpol and vpol  powers 




  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {
    TGraph * which_power = &powers[pol]; 
    for (int i = 0; i < NUM_SEAVEYS; i++) 
    {
      const TGraph * power = event->getFilteredGraph(i,(AnitaPol::AnitaPol_t)pol)->power(); 


      if (i == 0)
      {
        which_power->Set(power->GetN()); 
        memcpy(which_power->GetX(), power->GetX(), power->GetN() * sizeof(double)); 
        memset(which_power->GetY(), 0, power->GetN() * sizeof(double)); 
      }

      for (int j = 0; j < which_power->GetN(); j++) 
      {
          which_power->GetY()[j] += FFTtools::evalEvenGraph(power, which_power->GetX()[j])  / NUM_SEAVEYS ; 
      }
    }

     dBize(which_power); 
  }

  //compute means

  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {

    TGraph * baseline = pol == AnitaPol::kHorizontal ? hpol_avg : vpol_avg; 
    TGraph * power = &powers[pol]; 


    // tilt correction will make these change all the time 
    double baseline_mean= getMeanBetween(baseline,fmin,fmax); 
    double baseline_mean_first_half= getMeanBetween(baseline,fmin,(fmin+fmax)/2); 
    double baseline_mean_second_half= getMeanBetween(baseline,(fmin+fmax)/2, fmax); 
   
    double mean = getMeanBetween(power, fmin, fmax); 
    double mean_first_half = getMeanBetween(power, fmin, (fmin + fmax)/2); 
    double mean_second_half = getMeanBetween(power, (fmin + fmax)/2, fmax); 

    double delta_mean = mean - baseline_mean; 
    double delta_mean_first_half = mean_first_half - baseline_mean_first_half; 
    double delta_mean_second_half = mean_second_half - baseline_mean_second_half; 
    double slope = (delta_mean_first_half- delta_mean_second_half) / ( (fmax - fmin)/2); 

    //apply gain correction to power 
    for (int i = 0; i < power->GetN(); i++) 
    {
      power->GetY()[i] -= delta_mean; 
    }

    //apply tilt correction to baseline
    for (int i = getBin(fmin,baseline); i <= getBin(fmax,baseline); i++) 
    {
      baseline->GetY()[i] -= slope * (baseline->GetX()[i] - (fmax + fmin)/2 ); 
    }
  

    //find mean frequency ... is this really the mean?  seems more like the median to me... 

    double total_power = 0; 
    int power_min_bin = getBin(fmin,power); 
    int power_max_bin = getBin(fmax,power); 

    for (int i = power_min_bin; i <= power_max_bin; i++) 
    {
      total_power += TMath::Power(10, power->GetY()[i]/10); 
    }

    double cml_power = 0; //
    for (int i = power_min_bin; i <= power_max_bin; i++) 
    {
      if (cml_power >=  total_power/2)
      {
        mean_freq[pol] = power->GetX()[i]; 
        break; 
      }
      cml_power += TMath::Power(10, power->GetY()[i]/10); 
    }


    // see if any peaks are ndB above the baseline

    TGraph deltaMag(power->GetN(), power->GetX(), power->GetY()); 
    for (int i = power_min_bin; i <= power_max_bin; i++) 
    {
      deltaMag.GetY()[i]-= FFTtools::evalEvenGraph(baseline,power->GetX()[i]); 
    }

    strongest_cw[pol] = 0; 
    int max_freq_index = 0; 
    for (int i = 0; i < nfreq; i++) 
    {
       int max_index = TMath::MaxElement(power_max_bin- power_min_bin + 1,   deltaMag.GetY() + power_min_bin); 
       max_index += power_min_bin; 
       double val = deltaMag.GetY()[max_index]; 
       double f = deltaMag.GetX()[max_index]; 
//       printf("%f: %f\n", f, val); 

       if (val > dbCut)
       {
         freqs[pol][i] = f; 
         if (i == 0) 
         {
           strongest_cw[pol] = val; 
         }
       }
       else
       {
         max_freq_index = i; 
         break; 
       }

       for (int i = getBin(f-bw, &deltaMag); i <= getBin(f+bw, &deltaMag); i++) 
       {
         deltaMag.GetY()[i] = -1000; 
       }


    }

    for (int i =  max_freq_index; i < nfreq; i++) 
    {
      freqs[pol][i] = -1; 
    }

  }

  //now that we have the frequencies to cut, we should cut them!
  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {
    for (int freq = 0; freq < nfreq; freq++) 
    {
      if (freqs[pol][freq] == -1) break; 

      ComplicatedNotchFilter complicated(freqs[pol][freq]-bw, freqs[pol][freq]+bw, temperature, gain); 
//      printf("Filtering out %f\n", freqs[pol][freq]); 
      for (int i= 0; i < NUM_SEAVEYS; i++) 
      {
          complicated.processOne(getWf(event,i, AnitaPol::AnitaPol_t(pol))); 
      }
    }
  }
  
  
}



//
///// SINE SUBTRACT //// 
//


UCorrelator::SineSubtractFilter::SineSubtractFilter(double min_power_ratio, int max_failed_iter, double oversample_factor, int nfreq_bands, const double * fmin, const double * fmax, int nstored_freqs)
  : nstored_freqs(nstored_freqs) 
{
 desc_string.Form("SineSubtract with min_power_ratio=%f, max_failed_iter =%d, oversample_factor = %f and %d frequency bands", min_power_ratio, max_failed_iter,oversample_factor, nfreq_bands);
 if (nfreq_bands)
 {
   desc_string += ". The bands are: "; 
   for (int i = 0; i < nfreq_bands; i++)
   {
     desc_string += TString::Format ("(%f -> %f)", fmin[i], fmax[i]); 
   }
 }
 

 output_names.push_back(TString("total_power_removed")); 

 const char * polstr[2] = {"hpol","vpol"}; 
 for (int pol = 0; pol < 2; pol++)
 {
   for (int i = 0; i < NUM_SEAVEYS; i++) 
   {
     subs[pol][i] = new FFTtools::SineSubtract(max_failed_iter, min_power_ratio); 
     subs[pol][i]->setOversampleFactor(oversample_factor); 
     subs[pol][i]->setFreqLimits(nfreq_bands, fmin, fmax); 

   }

   output_names.push_back(TString::Format("niter_%s",polstr[pol])); 
   output_names.push_back(TString::Format("power_removed_%s",polstr[pol])); 

   for (int j = 0; j < nstored_freqs; j++)
   {
       output_names.push_back(TString::Format("freq%d_%s",j,polstr[pol])); 
       output_names.push_back(TString::Format("A%d_%s",j,polstr[pol])); 
       output_names.push_back(TString::Format("phase%d_%s",j,polstr[pol])); 
   }
 }
}


void UCorrelator::SineSubtractFilter::fillOutput(unsigned ui, double * vars) const 
{
//  printf("fillOutput(%u)\n",ui); 

  int i = int(ui);  //cheat; 
  if (i == 0)
  {
    double total_power_initial = 0; 
    double total_power_final = 0; 

    for (int pol = 0; pol < 2; pol++)
    {
      for (int i = 0; i < NUM_SEAVEYS; i++) 
      {
        const FFTtools::SineSubtractResult *r =  subs[pol][i]->getResult(); 
        double pinit = r->powers[0]; 
        double pfinal = r->powers[r->powers.size()-1]; 
        total_power_initial += pinit;
        total_power_final += pfinal; 
      }
    }
    *vars = 1. - total_power_final / total_power_initial; 
    return; 
  }

  int nfor_each_pol = (nOutputs()-1)/2; 

  int pol =  (i-1) < nfor_each_pol ? 0 : 1; 

  for (int j = 0; j < NUM_SEAVEYS; j++) 
  {
     const FFTtools::SineSubtractResult *r =  subs[pol][j]->getResult(); 


     if (i == 1 || i == 1 + nfor_each_pol) 
     {
       vars[j] = (double) r->freqs.size(); 
     }
     else if (i == 2 || i == 2 + nfor_each_pol) 
     {
       vars[j] = (double) r->powers[r->powers.size()-1] - r->powers[0]; 
     }
     else 
     {
       int ii =i -  3  -  (pol) * nfor_each_pol; 

       int freqi = ii /3; 
//       printf("%d\n", freqi); 

       if (freqi >= r->freqs.size())
       {
         vars[j] = -1; 
         continue; 
       }

       if (ii%3 ==0 ) 
       {
         vars[j] = r->freqs[freqi]; 
       }
       else if (ii %3 == 1) 
       {
         vars[j] = r->amps[0][freqi]; 
       }
       else
       {
         vars[j] = r->phases[0][freqi]; 
       }
     }
  }

}

void UCorrelator::SineSubtractFilter::process(FilteredAnitaEvent * ev) 
{

#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NUM_SEAVEYS; i++) 
  {
    for (int pol = 0; pol < 2; pol++)
    {
      AnalysisWaveform * wf = getWf(ev, i, AnitaPol::AnitaPol_t(pol)); 
      TGraph * g = wf->updateUneven(); 
      subs[pol][i]->subtractCW(1,&g, 1/2.6); 
    }
  }
}

UCorrelator::SineSubtractFilter::~SineSubtractFilter()
{
 for (int pol = 0; pol < 2; pol++)
 {
   for (int i = 0; i < NUM_SEAVEYS; i++) 
   {
     delete subs[pol][i]; 
   }
 }
}


