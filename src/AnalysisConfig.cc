#include "AnalysisConfig.h" 
#include "SystemResponse.h" 
#include "TFile.h" 
#include "TH2.h" 
#include "AnitaVersion.h"

static const char * peakfinders[] = {"Abby","Bicubic","Gaussian","QuadraticFit9","QuadraticFit16","QuadraticFit25", "Histogram" }; 
static const char * responses[] = {"None","SingleBRotter","IndividualBRotter","HarmSignalOnly", "TUFFs", "notches_260_375_0"}; 




#ifdef ENABLE_LIBCONFIG
#include "libconfig.h++"
#define LOOKUP(X) cfg.lookupValue(#X,X) 

void lookupEnum(libconfig::Config * cfg, const char * key, int * val, int N, const char ** allowed) 
{
  if (!cfg->exists(key)) return; 

  libconfig::Setting & set = cfg->lookup(key); 

  for (int i = 0; i < N; i++)
  {
    if (!strcasecmp(set, allowed[i]))
    {
      *val = i; 
      return; 
    }
  }

  fprintf(stderr,"Config Parse Error on line %d: %s must be one of: \n\t", set.getSourceLine(), key); 


  for (int i = 0; i < N; i++)
  {
    fprintf(stderr,"%s ", allowed[i]); 
  }

  fprintf(stderr, "\n But %s was given.", (const char * ) set); 

}


void UCorrelator::AnalysisConfig::loadFromFile(const char * config_file) 
{

  libconfig::Config cfg; 
  cfg.setAutoConvert(true); 
  cfg.readFile(config_file); 


  LOOKUP(correlator_nphi); 
  LOOKUP(correlator_ntheta); 
  LOOKUP(correlator_theta_lowest); 
  LOOKUP(correlator_theta_highest); 
  LOOKUP(enable_group_delay); 
  LOOKUP(zoomed_nphi); 
  LOOKUP(zoomed_ntheta); 
  LOOKUP(zoomed_dphi); 
  LOOKUP(zoomed_dtheta); 
  LOOKUP(zoomed_nant); 
  LOOKUP(combine_nantennas); 
  LOOKUP(combine_npad); 
  LOOKUP(saturation_threshold); 
  LOOKUP(peak_isolation_requirement); 
  LOOKUP(nmaxima); 
  LOOKUP(use_bin_center); 
  LOOKUP(bw_ndb);
  LOOKUP(spectral_fit_start);
  LOOKUP(spectral_fit_stop); 
  LOOKUP(noise_estimate_t0); 
  LOOKUP(noise_estimate_t1); 
  LOOKUP(response_npad); 
  LOOKUP(scale_by_cos_theta); 
  LOOKUP(use_offline_mask); 
  LOOKUP(max_mean_power_filtered); 
  LOOKUP(max_median_power_filtered); 
  LOOKUP(max_bottom_to_top_ratio); 
  LOOKUP(baseline_weight); 
  LOOKUP(max_peak_trigger_angle); 
  LOOKUP(min_peak_distance_from_unmasked); 
  LOOKUP(windowStokes);
  LOOKUP(stokesWindowLength);
  LOOKUP(set_bottom_first);
  LOOKUP(delay_to_center);
  LOOKUP(use_forced_trigger_rms); 
  LOOKUP(use_coherent_spectra); 
  LOOKUP(combine_t0); 
  LOOKUP(combine_t1); 
  const char * pols[] = {"horizontal", "vertical" }; 
  lookupEnum(&cfg, "start_pol", (int*) &start_pol, 2,pols); 
  lookupEnum(&cfg, "end_pol", (int*) &end_pol, 2,pols); 
  lookupEnum(&cfg, "fine_peak_finding_option", (int*) &fine_peak_finding_option, sizeof(peakfinders)/sizeof(char *), peakfinders); 
  lookupEnum(&cfg, "response_option", (int*) &response_option, sizeof(responses)/sizeof(char *), responses); 

}

#else
void UCorrelator::AnalysisConfig::loadFromFile(const char * config_file) 
{

  fprintf(stderr, "Not compiled with support for reading config files. You need libconfig for that.\n"); 
  (void) config_file; 
}
#endif


const int wais_hpol_time_offset[5] = {0,0,0,93,11000}; 
const int wais_vpol_time_offset[5] = {0,0,0,-99757,1000}; 
const int siple_hpol_time_offset = -41;
const int siple_vpol_time_offset = +328;

UCorrelator::AnalysisConfig::AnalysisConfig(const char * config) 
  : 
    wais_hpol(wais_hpol_time_offset[AnitaVersion::get()], 800e3, 1e3), 
    wais_vpol(wais_vpol_time_offset[AnitaVersion::get()], 800e3, 1e3), 
    siple_hpol(siple_hpol_time_offset, 800e3, 1e3), 
    siple_vpol(siple_vpol_time_offset, 800e3, 1e3)
{
  correlator_nphi = 180; 
  correlator_ntheta = 100; 
  correlator_theta_lowest = 60; 
  correlator_theta_highest = 40; 
  enable_group_delay = true; 
  use_offline_mask = true; 
  zoomed_nphi = 60; 
  zoomed_ntheta = 60; 
  zoomed_dphi = 0.5; 
  zoomed_dtheta = 0.5; 
  zoomed_nant = 12; 


  combine_nantennas = 12; 
  combine_npad = 3; 

  saturation_threshold = 2500; 

  start_pol = AnitaPol::kHorizontal; 
  end_pol = AnitaPol::kVertical; 

  peak_isolation_requirement = 10; 

  fine_peak_finding_option = FinePeakFindingQuadraticFit25; 
  nmaxima = 2; 
  use_bin_center = false; 
  scale_by_cos_theta = false; 
  max_peak_trigger_angle = 0; 
  min_peak_distance_from_unmasked = -1; 

  bw_ndb = 6; //let the bandwidth search stop immediately
  spectral_fit_start = 0.22; 
  spectral_fit_stop =1.2; 

  noise_estimate_t0 = 70; 
  noise_estimate_t1 = 100; 

  response_option = ResponseNone; 
  response_npad = 50; 

  max_mean_power_filtered = 1e6; 
  max_median_power_filtered = 1e6; 
  max_bottom_to_top_ratio = 5; 

  baseline_weight = 0; 

  if (config) loadFromFile(config);

  deconvolution_method = &AnitaResponse::kDefaultDeconvolution; 
  
  the_ldb_hist = 0; 

  ldb_max_run = 160; 

  windowStokes = true;
  stokesWindowLength = -1;
  
  set_bottom_first = false;

  delay_to_center = false;
  combine_t0 = -25; 
  combine_t1 = 125; 

}

static int nag = 0; 
TH2 * UCorrelator::AnalysisConfig::ldb_hist() const 
{
  if (the_ldb_hist) return the_ldb_hist; 

  TString fname; 
  fname.Form("%s/share/UCorrelator/ldbSelection.root", getenv("ANITA_UTIL_INSTALL_DIR")); 
  TFile f(fname); 

  if (!f.IsOpen())
  {
    if (nag++ < 3) 
      fprintf(stderr,"WARNING: data/ldbSelection.root doesn't exist. This is probably because you haven't generated it.\n. If you don't care about LDB pulsers, ignore this. Otherwise, you can generate one using the makeLDBSelection.C macro.\n"); 

    the_ldb_hist = 0; 
  }
  else
  {
    the_ldb_hist = (TH2*) f.Get("ldbHist")->Clone("LDBHistogram"); 
    the_ldb_hist->SetDirectory(0); 
  }
  return the_ldb_hist; 
}


const char * UCorrelator::AnalysisConfig::getPeakFindingString(FinePeakFindingOption_t opt) 
{
  return peakfinders[opt]; 
}

const char * UCorrelator::AnalysisConfig::getResponseString(ResponseOption_t opt) 
{
  return responses[opt]; 
}


UCorrelator::AnalysisConfig::~AnalysisConfig()
{
  if (the_ldb_hist) delete the_ldb_hist; 
}

