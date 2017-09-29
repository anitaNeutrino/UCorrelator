#include "FFTtools.h" 
#include "FilterStrategy.h"
#include "UCFilters.h" 
#include "SpectrumAverage.h" 
#include "GeomFilter.h" 
#include "BasicFilters.h" 
#include "ResponseManager.h" 
#include "SystemResponse.h" 
#include "TGraph.h" 
#include "AnalysisWaveform.h" 
#include "RawAnitaHeader.h"
#include "TRandom.h"
#include "TaperFilter.h" 
#include "Adu5Pat.h"
#include "DigitalFilter.h"
#include <math.h>// for isnan
#include "AntennaPositions.h"
#include "SineSubtractCache.h"

#ifdef UCORRELATOR_OPENMP
#include "omp.h"
#endif




FilterStrategy * UCorrelator::getStrategyWithKey (const char * key) 
{
  FilterStrategy * s = new FilterStrategy; 
  fillStrategyWithKey(s,key); 
  return s; 
}


static std::map<const char *, TString> key_descs; 
static UCorrelator::SpectrumAverageLoader avgldr; 
static AnitaResponse::ResponseManager * responseManager = 0; 
static AnitaResponse::AllPassDeconvolution allpass; 
// static UCorrelator::SpectrumAverageLoader avgldr; 
// static UCorrelator::ResponseManager * responseManager = 0; 
// static UCorrelator::AllPassDeconvolution allpass; 


static TString desc; 
const char * UCorrelator::fillStrategyWithKey(FilterStrategy * fillme, const char * key) 
{
  TString tokens(key); 
  TString tok; 
  Ssiz_t from = 0; 

  desc = ""; 
  bool need_description = key_descs.count(key); 

  while (tokens.Tokenize(tok, from,"+"))
  {

    if (desc.Length() > 0) desc += " + ";  

    if (strcasestr(tok.Data(),"sinsub_") == tok.Data())
    {
      int mpr; int iter; 
      if(sscanf(tok.Data(),"sinsub_%02d_%d",&mpr,&iter)!=2) 
      {
        fprintf(stderr,"Problem with token: %s\n", tok.Data()); 
      }
      else
      {
        SineSubtractFilter * ssf = new SineSubtractFilter(mpr/100.,iter);

        /** now check for additional options */ 
        if (const char * substr =strstr(tok.Data(),"_ad"))
        {
          int exp; 
          if (sscanf(substr,"_ad_%d", &exp) == 1)
          {
            double expf = exp; 
            while (expf > 10) expf /=10; 
            ssf->makeAdaptive(&avgldr,expf); 
          }
          else
          {
            fprintf(stderr, "Problem with _ad option in token: %s\n", tok.Data()); 
          }
        }

        if (const char * substr = strstr(tok.Data(),"_env"))
        {
          const char * parstr = 0; 
          FFTtools::SineSubtract::EnvelopeOption env = FFTtools::SineSubtract::ENV_NONE; 
          double pars[3]; 
          int nargs =0;
          if (strstr(substr,"_env_peak")==substr)
          {
            parstr = substr + strlen("_env_peak"); 
            env = FFTtools::SineSubtract::ENV_PEAK;  
          }
          else if (strstr(substr,"_env_rms")==substr)
          {
            parstr = substr + strlen("_env_rms"); 
            env = FFTtools::SineSubtract::ENV_RMS;  
          }
          else if (strstr(substr,"_env_hilbert")==substr)
          {
            parstr = substr + strlen("_env_hilbert"); 
            env = FFTtools::SineSubtract::ENV_RMS;  
          }

          /* try to parse arguments */ 

          if (parstr) 
          {
            nargs = sscanf(parstr,"_%lg_%lg_%lg",pars,pars+1,pars+2); 
          }
          ssf->setEnvelopeOption(env, nargs ? pars : 0);
        }

        if (const char * substr = strstr(tok.Data(),"_peakfind"))
        {
          const char * parstr = 0; 
          FFTtools::SineSubtract::PeakFindingOption peak = FFTtools::SineSubtract::NEIGHBORFACTOR; 
          double pars[4]; 
          int nargs = 0;
          if (strstr(substr,"_peakfind_global")==substr)
          {
            parstr = substr + strlen("_peakfind_global"); 
            peak = FFTtools::SineSubtract::GLOBALMAX;  
          }
          else if (strstr(substr,"_peakfind_neighbor")==substr)
          {
            parstr = substr + strlen("_peakfind_neighbor"); 
            peak = FFTtools::SineSubtract::NEIGHBORFACTOR;  
          }
          else if (strstr(substr,"_peakfind_tspectrum")==substr)
          {
            parstr = substr + strlen("_peakfind_tspectrum"); 
            peak = FFTtools::SineSubtract::TSPECTRUM;  
          }
          else if (strstr(substr,"_peakfind_savgolsub")==substr)
          {
            parstr = substr + strlen("_peakfind_savgolsub"); 
            peak = FFTtools::SineSubtract::SAVGOLSUB;  
          }

          /* try to parse arguments */ 

          if (parstr) 
          {
            nargs = sscanf(parstr,"_%lg_%lg_%lg_%lg",pars,pars+1,pars+2, pars+3); 
          }
          ssf->setPeakFindingOption(peak, nargs ? pars : 0);
        }
          


        if (fillme) 
        {
          fillme->addOperation(ssf); 
          // fillme->addOperation(ssf,true);// true for having output. 
        }

        if (need_description)
        {
          desc+= ssf->description(); 
        }

        if (!fillme)
        {
          delete ssf; 
        }
      }
    }
    /** shortcutoption */ 
    else if (strcasestr(tok.Data(),"adsinsub_"))
    {
      int exp; int mpr; int iter; 
      if(sscanf(tok.Data(),"adsinsub_%d_%02d_%d",&exp,&mpr,&iter)!=3) 
      {
        fprintf(stderr,"Problem with token: %s\n", tok.Data()); 
      }
      else
      {

        double expf = exp; 
        while (expf > 10) expf /=10; 

        UCorrelator::SineSubtractFilter * ssf = new UCorrelator::SineSubtractFilter(mpr/100.,iter); 
        if (fillme) 
        {
          ssf->makeAdaptive(&avgldr,expf); 
          fillme->addOperation(ssf); 
        }
        if (need_description)
        {
          desc+= ssf->description(); 
        }

        if (!fillme)
          delete ssf;
      }
    }

    else if (strcasestr(tok.Data(),"butter_")) 
    {
      int thresh, order; 
      if(sscanf(tok.Data(),"butter_%d_%d",&thresh,&order)!=2) 
      {
        fprintf(stderr,"Problem with token: %s\n", tok.Data()); 
      }
      else
      {
        double threshf = thresh; 
        while (threshf > 10) threshf/=10; 

        if (fillme) 
        {
          fillme->addOperation(new AdaptiveButterworthFilter(&avgldr,threshf,order)); 
        }
        if (need_description) 
        {
          desc += TString::Format("(Adaptive Butterworth filter with order %d, peakiness threshold %g)",order,threshf);
        }
      }
    }
    else if (strcasestr(tok.Data(),"minphase_")) 
    {
      int exp;
      if(sscanf(tok.Data(),"minphase_%d",&exp)!=1) 
      {
        fprintf(stderr,"Problem with token: %s\n", tok.Data()); 
      }
      else
      {
        double expf = exp; 
        while (expf > 10 )expf/=10; 

        if (fillme) 
          fillme->addOperation(new AdaptiveMinimumPhaseFilter(&avgldr,-expf)); 

        if (need_description)
        {
          desc += TString::Format("(Adaptive minimum phase filter with exp %g)", expf); 
        }
      }
    }
    else if (strcasestr(tok.Data(),"brickwall_")) 
    {
      int thresh;int fill;
      if(sscanf(tok.Data(),"brickwall_%d_%d",&thresh,&fill)!=2) 
      {
        fprintf(stderr,"Problem with token: %s\n", tok.Data()); 
      }
      else
      {
        double threshf = thresh; 
        while (threshf > 10 )threshf/=10; 

        if (fillme) 
          fillme->addOperation(new AdaptiveBrickWallFilter(&avgldr,threshf,fill)); 

        if (need_description)
        {
          desc += TString::Format("(Adaptive brick wall filter with thresh %g, filNotch=%d)", threshf,fill); 
        }
      }
    }


    else if (strcasestr(tok.Data(),"deconv")) 
    {
      //TODO: A4 support 
      if (!responseManager && fillme ) 
      {
        responseManager = new AnitaResponse::ResponseManager("IndividualBRotter",3); 
      }

      if (need_description) 
      {
        desc += TString::Format("(deconv)"); 
      }

      if (fillme) 
        fillme->addOperation( new AnitaResponse::DeconvolveFilter(responseManager, &allpass)); 
    }

    else if (strcasestr(tok.Data(),"gaustaper_")) 
    {

      double mean,nsigma; 


      if(sscanf(tok.Data(),"gaustaper_%lg_%lg",&mean,&nsigma)!=2) 
      {
        fprintf(stderr,"Problem with token: %s\n", tok.Data()); 
      }
      else
      {
        if (need_description) 
        {
          desc += TString::Format("(Gaussian taper, distance=%lg, nsigma=%lg)",mean,nsigma); 
        }


        if (fillme) 
          fillme->addOperation(new GaussianTaper(mean,nsigma)); 

      }
 
    }


    else if (strcasestr(tok.Data(),"geom"))
    {
      if (fillme) 
      {
        std::vector<std::vector<TGraphAligned* > > noise(2, std::vector<TGraphAligned*>(48)); 
        for (int ant = 0; ant < 48; ant++)
        {
          for (int pol = 0; pol < 2; pol++) 
          {
            TH1 * avg = UCorrelator::SpectrumAverage::defaultThermal()->getSpectrumPercentile(AnitaPol::AnitaPol_t(pol), ant,0.5,true); 
            noise[pol][ant] = new TGraphAligned(avg->GetNbinsX()); 
            for (int i = 0; i < avg->GetNbinsX(); i++) 
            {
              noise[pol][ant]->GetX()[i] = avg->GetBinLowEdge(i+1); 
              noise[pol][ant]->GetY()[i] = avg->GetBinContent(i+1); 
            }
            delete avg; 
          }
        }
        fillme->addOperation(new GeometricFilter(noise)); 
      }
      if (need_description)
      {
        desc += "(GeometricFilter)"; 
      }
    }


    else if (strcasestr(tok.Data(),"abby"))
    {
      if (fillme) 
      {
        fillme->addOperation(new SimplePassBandFilter(0.2, 1.2)); 


        fillme->addOperation(new ConditionalFilterOperation
                                     (  new ComplicatedNotchFilter(0.26-0.026, 0.26 + 0.026), 
                                            &UCorrelator::antennaIsNorthFacing, "_north", "if is facing north",true
                                     )
                             );  

        fillme->addOperation(new AdaptiveFilterAbby(2, getenv("UCORRELATOR_BASELINE_DIR"))); 

      }
      if (need_description)
      {
        desc += "(AbbyFilter)"; 
      }
    }
  }
  


  if (AnitaVersion::get()==3) 
  {
    //TODO: make this more configurable
    if (fillme) 
    {
      fillme->addOperation(new ALFAFilter); 
    }
    if (need_description) 
    {
      desc += "(ALFA Filter) "; 
    }
  }

  if (need_description) 
    key_descs[key] = desc; 
  else
    desc = key_descs[key]; 


  return desc.Data(); 
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

const char * UCorrelator::AdaptiveFilterAbby::outputName(unsigned i) const
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


void UCorrelator::AdaptiveFilterAbby::fillOutput(unsigned  i, double *vals) const
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

void UCorrelator::AdaptiveFilterAbby::process(FilteredAnitaEvent * event) 
{
//  printf("AdaptiveFilterAbby::process!\n"); 
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



static bool use_sine_sub_cache = false;
void UCorrelator::SineSubtractFilter::setUseCache(bool uc){
  use_sine_sub_cache = uc;
}



UCorrelator::SineSubtractFilter::SineSubtractFilter(double min_power_ratio, int max_failed_iter,  int nfreq_bands, const double * fmin, const double * fmax, int nstored_freqs)
    : min_power_ratio(min_power_ratio), spec(0), last_t(0), nstored_freqs(nstored_freqs), adaptive_exp(1), max_failed(max_failed_iter), sine_sub_cache(NULL) 
{

 memset(reduction,0,sizeof(reduction)); 

 desc_string.Form("SineSubtract with min_power_ratio=%f, max_failed_iter =%d, and %d frequency bands", min_power_ratio, max_failed_iter, nfreq_bands);
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

   for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
     AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
     for(int ant=0; ant < NUM_SEAVEYS; ant++){
       cached_ssr[pol][ant] = NULL;
     }
   }
}


void UCorrelator::SineSubtractFilter::setPeakFindingOption(FFTtools::SineSubtract::PeakFindingOption option, const double * params)
{
 for (int pol = 0; pol < 2; pol++)
 {
   for (int i = 0; i < NUM_SEAVEYS; i++) 
   {
     subs[pol][i]->setPeakFindingOption(option,params); 
   }
 }
}

void UCorrelator::SineSubtractFilter::setEnvelopeOption(FFTtools::SineSubtract::EnvelopeOption option, const double * params)
{
 for (int pol = 0; pol < 2; pol++)
 {
   for (int i = 0; i < NUM_SEAVEYS; i++) 
   {
     subs[pol][i]->setEnvelopeOption(option,params); 
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

       if (freqi >= (int) r->freqs.size())
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


void UCorrelator::SineSubtractFilter::refresh_cache(UInt_t eventNumber){
  if(use_sine_sub_cache){
    if(!sine_sub_cache){
      sine_sub_cache = new SineSubtractCache(this->description());
    }
  }
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      if(use_sine_sub_cache){
        cached_ssr[pol][ant] = sine_sub_cache->getResult(eventNumber, pol, ant);
      }
      else{
        cached_ssr[pol][ant] = NULL;
      }
    }
  }
}

void UCorrelator::SineSubtractFilter::process(FilteredAnitaEvent * ev) 
{
 const RawAnitaHeader * h = ev->getHeader(); 
 refresh_cache(h->eventNumber);

#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < 2*NUM_SEAVEYS; i++) 
  {
      int pol = i % 2; 
      AnalysisWaveform * wf = getWf(ev, i/2, AnitaPol::AnitaPol_t(pol)); 
      processOne(wf,h,i/2,pol);
  }

  last_t = ev->getHeader()->triggerTime; //ev->getHeader()->triggerTime; 
}

void UCorrelator::SineSubtractFilter::processOne(AnalysisWaveform *wf, const RawAnitaHeader * header, int i, int pol)
{

  if (spec) 
  {
    // use peakiness to tune aggressivness
    if (header->triggerTime > last_t)
    {
      TH2 * peaky = (TH2*) spec->avg(header->triggerTime)->getPeakiness(AnitaPol::AnitaPol_t(pol), i); 

      if (reduction[pol][i])
        delete reduction[pol][i]; 

      reduction[pol][i] = new TGraph(peaky->GetNbinsX()); 

      for (int j = 0; j < reduction[pol][i]->GetN(); j++) 
      {
        reduction[pol][i]->GetX()[j] = peaky->GetXaxis()->GetBinLowEdge(j+1); 
        double t = header->triggerTime + header->triggerTimeNs * 1e-9;

        if (header->triggerTime < peaky->GetYaxis()->GetXmin())
        {
             fprintf(stderr,"Warning, time %g before first point in histogram: %g\n", t, peaky->GetYaxis()->GetXmin()); 
             t = peaky->GetYaxis()->GetBinCenter(1); 
        }

        if (header->triggerTime > peaky->GetYaxis()->GetXmax())
        {
             fprintf(stderr,"Warning, time %g after last point in histogram: %g\n", t, peaky->GetYaxis()->GetXmax()); 
             t = peaky->GetYaxis()->GetBinCenter(peaky->GetNbinsY()); 
        }


        double how_peaky = peaky->Interpolate(peaky->GetXaxis()->GetBinCenter(j+1), t); 
        if (how_peaky < 1) how_peaky = 1; 
        reduction[pol][i]->GetY()[j] = min_power_ratio/TMath::Power(how_peaky,adaptive_exp); 
      }
    }

    subs[pol][i]->setMinPowerReduction(reduction[pol][i]); 
  }
  else
  {
    subs[pol][i]->setMinPowerReduction(min_power_ratio); 
  }


  TGraph * g = wf->updateUneven();
  if (g->GetRMS(2) > 0){
    // subs[pol][i]->subtractCW(1,&g, 1/2.6, NULL, cached_ssr[pol][i]);

    // if(cached_ssr[pol][i]){
    //   TGraph g2 = *g;
    //   TGraph* gr2p = &g2;
    //   std::cout << pol << "\t" << i << "\t" << cached_ssr[pol][i]->powers.at(0) << "\t";
    //   subs[pol][i]->subtractCW(1, &gr2p, 1/2.6, NULL, cached_ssr[pol][i]);
    //   subs[pol][i]->subtractCW(1, &g, 1/2.6, NULL);
    //   auto r = subs[pol][i]->getResult();
    //   std::cerr << r->powers.at(0) << std::endl;
    // }
    // else{
    subs[pol][i]->subtractCW(1,&g, 1/2.6, NULL, cached_ssr[pol][i]);
    // }
  }
}

UCorrelator::SineSubtractFilter::~SineSubtractFilter()
{
 for (int pol = 0; pol < 2; pol++)
 {
   for (int i = 0; i < NUM_SEAVEYS; i++) 
   {
     delete subs[pol][i]; 

     if (reduction[pol][i]) 
       delete reduction[pol][i]; 
   }
 }
}


void UCorrelator::SineSubtractFilter::setVerbose(bool set) 
{
 for (int pol = 0; pol < 2; pol++)
 {
   for (int i = 0; i < NUM_SEAVEYS; i++) 
   {
     subs[pol][i]->setVerbose(set); 
   }
 }


}


void UCorrelator::SineSubtractFilter::setInteractive(bool set) 
{
 for (int pol = 0; pol < 2; pol++)
 {
   for (int i = 0; i < NUM_SEAVEYS; i++) 
   {
     subs[pol][i]->setStore(set); 
   }
 }


}


void UCorrelator::SineSubtractFilter::makeAdaptive(const SpectrumAverageLoader * s, double peak_exp) 
{
  spec = s; 
  adaptive_exp = peak_exp; 
  //TODO: include frequency bands
  desc_string.Form("Adaptive SineSubtract with adaptive exponent=%f, min_power_ratio=%f, max_failed_iter =%d", peak_exp, min_power_ratio, max_failed);
}


//Combined Sine Subtract 

UCorrelator::CombinedSineSubtractFilter::CombinedSineSubtractFilter(double min_power_ratio, int max_failed_iter, int nfreq_bands, const double * fmin, const double * fmax, int nstored_freqs)
  : nstored_freqs(nstored_freqs) 
{
 desc_string.Form("CombinedSineSubtract with min_power_ratio=%f, max_failed_iter =%d, and %d frequency bands", min_power_ratio, max_failed_iter, nfreq_bands);
 if (nfreq_bands)
 {
   desc_string += ". The bands are: "; 
   for (int i = 0; i < nfreq_bands; i++)
   {
     desc_string += TString::Format ("(%f -> %f)", fmin[i], fmax[i]); 
   }
 }
 

 output_names.push_back(TString("total_power_removed")); 

 for (int i = 0; i < NUM_PHI; i++) 
 {
   for (int ipol = 0; ipol < 2; ipol++)
   {
     subs[i][ipol] = new FFTtools::SineSubtract(max_failed_iter, min_power_ratio); 
     subs[i][ipol]->setFreqLimits(nfreq_bands, fmin, fmax); 
   }
 }

 output_names.push_back("niter"); 
 output_names.push_back("power_removed"); 

 for (int j = 0; j < nstored_freqs; j++)
 {
     output_names.push_back(TString::Format("freq%d",j)); 
     output_names.push_back(TString::Format("A%d",j)); 
     output_names.push_back(TString::Format("phase%d",j)); 
 }
}


void UCorrelator::CombinedSineSubtractFilter::fillOutput(unsigned ui, double * vars) const 
{
//  printf("fillOutput(%u)\n",ui); 

  int i = int(ui);  //cheat; 
  if (i == 0)
  {
    double total_power_initial = 0; 
    double total_power_final = 0; 

    for (int i = 0; i < NUM_PHI; i++) 
    {
      for (int ipol = 0; ipol < 2; ipol++) 
      {
        const FFTtools::SineSubtractResult *r =  subs[i][ipol]->getResult(); 
        double pinit = r->powers[0]; 
        double pfinal = r->powers[r->powers.size()-1]; 
        total_power_initial += pinit;
        total_power_final += pfinal; 
      }
    }

    *vars = 1. - total_power_final / total_power_initial; 
    return; 
  }


  for (int j = 0; j < NUM_PHI; j++) 
  {
    for (int ipol = 0; ipol < 2; ipol++) 
    {

      const FFTtools::SineSubtractResult *r =  subs[j][ipol]->getResult(); 

      if (i == 1) 
      {
        vars[j + ipol * NUM_PHI] = (double) r->freqs.size(); 
      }
      else if (i == 2) 
      {
        vars[j + ipol * NUM_PHI] = (double) r->powers[r->powers.size()-1] - r->powers[0]; 
      }
      else 
      {
        int ii =i -  3; 

        int freqi = ii /3; 
//        printf("%d\n", freqi); 

        if (freqi >= (int) r->freqs.size())
        {
          vars[j] = -1; 
          continue; 
        }

        if (ii%3 ==0 ) 
        {
          vars[j + ipol * NUM_PHI] = r->freqs[freqi]; 
        }

        else 
        {
          for (unsigned k = 0; k < 3; k++)
          {
            int jj = j *3 +  k + ipol * NUM_PHI * 3; 

            if (ii %3 == 1) 
            {
              vars[jj] = r->amps.size() <= k ? 0 : r->amps[k][freqi]; 
            }
            else
            {
              vars[jj] = r->phases.size() <=k ? 0 : r->phases[k][freqi]; 
            }
          }
        }
      }
    }
  }
}

unsigned UCorrelator::CombinedSineSubtractFilter::outputLength(unsigned i) const
{

  if ( i == 0) return 1;  //total power removed
  if (i == 1 || i == 2) return NUM_PHI*2;  //power removed and niter per phi sector/ pol

  i-=3; 
  if (i % 3 == 0) return NUM_PHI*2; //frequency per phi sector /pol

  return NUM_PHI * 6; //phase and amp per antenna 
}

void UCorrelator::CombinedSineSubtractFilter::process(FilteredAnitaEvent * ev) 
{

#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NUM_PHI; i++) 
  {
      for (int pol = 0; pol < 2; pol++)
      {

        TGraph *gs[3]; //nrings * npol 
        int ng = 0; 

        for (int ring = 0; ring < 3; ring++)
        {
          AnalysisWaveform * wf = getWf(ev, i+ring*NUM_PHI, AnitaPol::AnitaPol_t(pol)); 
          TGraph * g = wf->updateUneven(); 
          if (g->GetRMS(2) > 0) 
            gs[ng++] = g; 
        }
        subs[i][pol]->subtractCW(ng,gs, 1/2.6); //NULL, cached_ssr[pol][i]); 
      }
   }
}

UCorrelator::CombinedSineSubtractFilter::~CombinedSineSubtractFilter()
{
   for (int i = 0; i < NUM_PHI; i++) 
   {
     for (int ipol = 0; ipol < 2; ipol++) 
     {
       delete subs[i][ipol]; 
     }
   }
}


void UCorrelator::CombinedSineSubtractFilter::setInteractive(bool set) 
{
   for (int i = 0; i < NUM_PHI; i++) 
   {
     for (int ipol = 0; ipol < 2; ipol++) 
     {
       subs[i][ipol]->setStore(set); 
     }
   }
}








UCorrelator::AdaptiveMinimumPhaseFilter::AdaptiveMinimumPhaseFilter(const SpectrumAverageLoader *avg, double exponent,int npad)
 : avg(avg),npad(npad), exponent(exponent), last_bin(-1) 
{
  desc_string.Form("Adaptive Minimum Phase Filter with SpectrumAverageLoader(%d) and exponent %g",  avg->getNsecs(), exponent); 
  memset(filt,0,sizeof(filt)); 
  memset(size,0,sizeof(size)); 
}

UCorrelator::AdaptiveMinimumPhaseFilter::~AdaptiveMinimumPhaseFilter() 
{
  for (int i = 0; i < 2; i++) 
  {
    for (int j = 0; j < NUM_SEAVEYS; j++) 
    {
      if (filt[i][j]) delete filt[i][j]; 
    }
  }
}


void UCorrelator::AdaptiveMinimumPhaseFilter::process(FilteredAnitaEvent * ev) 
{

  double t = ev->getHeader()->triggerTime + ev->getHeader()->triggerTimeNs*1e-9; 
  int bin = avg->avg(t)->getPeakiness(AnitaPol::kHorizontal,0)->GetYaxis()->FindBin(t); 


  for (int ipol = 0; ipol < 2; ipol++) 
  {
    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(ipol); 
    for (int i = 0; i < NUM_SEAVEYS; i++) 
    {
      AnalysisWaveform * wf = getWf(ev,i,pol); 


      //have to recalculate the filter
      if (last_bin != bin || wf->Neven() != size[ipol][i]) 
      {

        if (filt[ipol][i]) delete filt[ipol][i]; 
        double tmid = avg->avg(t)->getPeakiness(pol,i)->GetYaxis()->GetBinCenter(bin); 
        size[ipol][i] = wf->Neven(); 
        int fft_size = (size[ipol][i] * (1+npad))/2 +1; 

        double G[fft_size]; 
        for (int j = 0; j < fft_size; j++)
        {
          double f = 0.01/(1.+npad) *j; //TODO
//          printf("%g %g %g\n",f,tmid,t); 
          double peaky = ((TH2*)avg->avg(t)->getPeakiness(pol,i))->Interpolate(f,tmid);  //why isn't this const? ?!?? 
          if (peaky < 1 || std::isnan(peaky) || f < 0.16) peaky = 1; 
          G[j] = TMath::Power(peaky, exponent)*( j > 0 && j < fft_size-1 ? sqrt(2) : 1); 
        }

        filt[ipol][i] = FFTtools::makeMinimumPhase(fft_size,G); 
      }

      //zero pad? 
      if (npad) 
        wf->padEven(npad); 

      FFTWComplex * update = wf->updateFreq(); 
      for (int j = 0; j <size[ipol][i]*(1+npad)/2+1; j++) update[j] *= filt[ipol][i][j]; 

      //truncate
      if (npad) 
        wf->forceEvenSize(wf->Neven()/(1.+npad)); 

    }
  }

  last_bin = bin; 

}

TGraph * UCorrelator::AdaptiveMinimumPhaseFilter::getCurrentFilterTimeDomain(AnitaPol::AnitaPol_t pol, int i) const 
{
  int  N = size[pol][i] * (1+npad); 
  double * y = FFTtools::doInvFFT(N,filt[pol][i]); 

  TGraph * g = new TGraph(N); 

  for (int j = 0; j < N; j++)
  {
    g->GetX()[j]=j/(2.6); 
    g->GetY()[j]=y[j]; 
  }
  delete y; 
  return g; 
}


TGraph * UCorrelator::AdaptiveMinimumPhaseFilter::getCurrentFilterPower(AnitaPol::AnitaPol_t pol, int i) const 
{
  //TODO make right size

  int N = size[pol][i] * (1+npad) / 2 + 1; 
  
  TGraph * g = new TGraph(N); 

  for (int j = 0; j < N; j++)
  {
    g->GetX()[j]=j * 0.01 / (1+npad); 
    g->GetY()[j]=filt[pol][i][j].getAbsSq(); 
  }
  return g; 
}


  ////////////////////////
  ////AdaptiveBrickWall/// 
  ////////////////////////


static int n_adaptive_brickwall = 0; 
UCorrelator::AdaptiveBrickWallFilter::AdaptiveBrickWallFilter(const SpectrumAverageLoader * spec, double thresh, bool fillNotch) 
  : avg(spec), threshold(thresh), fill(fillNotch) 
{

  desc_string.Form("AdaptiveBrickWallFilter with Threshold=%g, fillNotch=%d", thresh,fillNotch); 
  last_bin = -1; 
  memset(sp,0,sizeof(sp)); 
  instance = n_adaptive_brickwall++; 
}

UCorrelator::AdaptiveBrickWallFilter::~AdaptiveBrickWallFilter()
{

  for (int i = 0; i < 2; i++) 
  {
    for (int j = 0; j < NUM_SEAVEYS; j++) 
    {
      if (sp[i][j]) delete sp[i][j]; 
    }
  }
}

void UCorrelator::AdaptiveBrickWallFilter::process(FilteredAnitaEvent *ev)
{

  double t = ev->getHeader()->triggerTime + ev->getHeader()->triggerTimeNs*1e-9; 
//  double t= ev->getHeader()->realTime; //until icemc fixes this
  int bin = avg->avg(t)->getPeakiness(AnitaPol::kHorizontal,0)->GetYaxis()->FindBin(t); 



  if (last_bin != bin) /* This has to happen in the main thread since ProjectionX() grabs some Mutex */ 
  {
    for (int ipol = 0; ipol < 2; ipol++) 
    {
      AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(ipol); 

      for (int i = 0; i < NUM_SEAVEYS; i++) 
      {
          if (sp[ipol][i]) 
          {
            delete sp[ipol][i]; 
          }
          sp[ipol][i] = avg->avg(t)->getPeakiness(pol,i)->ProjectionX(TString::Format("sp_%d_%d_%d",ipol,i,instance), bin,bin);  
          sp[ipol][i]->SetDirectory(0); 
      }
    }
  }

  for (int ipol = 0; ipol < 2; ipol++) 
  {
      AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(ipol); 
#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for 
#endif
      for (int i = 0; i < NUM_SEAVEYS; i++) 
      {
        AnalysisWaveform * wf = getWf(ev,i,pol); 
        int nf =wf->Nfreq(); 
        double df = wf->deltaF(); 
        FFTWComplex * fft = wf->updateFreq(); 
        for (int j =0; j < nf; j++)
        {
          double f = j * df; 
          double peaky = sp[ipol][i]-> Interpolate(f); 

          if (peaky > threshold) 
          {

            if (fill)
            {
              //random phase
              double phase = gRandom->Uniform(0, 2*M_PI); 
              double mag = fft[j].getAbs()/peaky; 
              fft[j].re = mag * cos(phase); 
              fft[j].im = mag * sin(phase); 
            }
            else
            {
              fft[j].re = 0; 
              fft[j].im = 0; 
            }
          }
        }
      }
  }
     
  last_bin = bin; 

}



UCorrelator::AdaptiveButterworthFilter::AdaptiveButterworthFilter(const SpectrumAverageLoader *avg,
                                                                  double peakiness_threshold ,
                                                                  int order , double width) 
            : avg(avg), threshold(peakiness_threshold), order(order), width(width)  
{ 
  desc_string.Form("AdaptiveButterworthFilter with SpectrumAverageLoader(%d), th=%g, order=%d, width=%g)",
                    avg->getNsecs(), threshold, order, width);
}

void UCorrelator::AdaptiveButterworthFilter::process(FilteredAnitaEvent * ev) 
{

  double t = ev->getHeader()->triggerTime + ev->getHeader()->triggerTimeNs*1e-9; 
//  double t = ev->getHeader()->realTime; //until icemc fixes this 
  int bin = avg->avg(t)->getPeakiness(AnitaPol::kHorizontal,0)->GetYaxis()->FindBin(t); 

  for (int ipol = 0; ipol < 2; ipol++) 
  {
    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(ipol); 
    for (int i = 0; i < NUM_SEAVEYS; i++) 
    {

      if (last_bin != bin) 
      {
        filters[ipol][i].clear();
        TH1 * proj = avg->avg(t)->getPeakiness(pol,i)->ProjectionX("tmp", bin,bin);  

        for (int j =1; j <= proj->GetNbinsX(); j++) 
        {
          if (proj->GetBinContent(j) < threshold)
          {
            continue; 
          }

           int low = j; 
           while (proj->GetBinContent(j++) > threshold) continue; 
           int high = j; 

           double fmax = proj->GetXaxis()->GetXmax(); 
           double freq_low = proj->GetBinLowEdge(low); 
           double freq_high = proj->GetBinLowEdge(high); 
           double freq_center = (freq_low+freq_high)/2; 
           double freq_width = width + freq_high - freq_low; 
           filters[ipol][i].add(new FFTtools::ButterworthFilter(FFTtools::NOTCH, order, freq_center/fmax, freq_width/(2*fmax)),true); 
        }

        delete proj; 
      }



      AnalysisWaveform * wf = getWf(ev,i,pol); 
      TGraph * update = wf->updateEven(); 
      filters[pol][i].filterGraph(update); 
     
    }
  }

  last_bin = bin; 

}


