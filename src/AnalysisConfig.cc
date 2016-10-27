#include "AnalysisConfig.h" 
#include "SystemResponse.h" 

static const char * peakfinders[] = {"Abby","Bicubic","Gaussian","QuadraticFit9","QuadraticFit16","QuadraticFit25", "Histogram" }; 
static const char * responses[] = {"None","SingleBRotter","IndividualBRotter","HarmSignalOnly"}; 




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
  LOOKUP(combine_unfiltered); 
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

  const char * pols[] = {"horizontal", "vertical" }; 
  lookupEnum(&cfg, "start_pol", (int*) &start_pol, 2,pols); 
  lookupEnum(&cfg, "end_pol", (int*) &end_pol, 2,pols); 
  lookupEnum(&cfg, "fine_peak_finding_option", (int*) &fine_peak_finding_option, sizeof(peakfinders)/sizeof(char *), peakfinders); 
  lookupEnum(&cfg, "response_option", (int*) &response_option, sizeof(responses)/sizeof(char *), responses); 

}

#else
void AnalysisConfig::loadFromFile(const char * config_file) 
{

  fprintf(stderr, "Not compiled with support for reading config files. You need libconfig for that.\n"); 
  (void) config_file; 
}
#endif


const int wais_hpol_time_offset = 93; 
const int wais_vpol_time_offset = -99757; 
const int ldb_hpol_time_offset = 0;  //TODO
const int ldb_vpol_time_offset = 0; 

UCorrelator::AnalysisConfig::AnalysisConfig(const char * config) 
  : 
    wais_hpol(wais_hpol_time_offset, 800e3, 1e3), 
    wais_vpol(wais_vpol_time_offset, 800e3, 1e3), 
    ldb_hpol(ldb_hpol_time_offset, 800e3, 1e3), 
    ldb_vpol(ldb_vpol_time_offset, 800e3, 1e3) 
{
  correlator_nphi = 180; 
  correlator_ntheta = 90; 
  correlator_theta_lowest = 60; 
  correlator_theta_highest = 25; 
  enable_group_delay = true; 


  zoomed_nphi = 60; 
  zoomed_ntheta = 60; 
  zoomed_dphi = 0.5; 
  zoomed_dtheta = 0.5; 
  zoomed_nant = 0; 


  combine_nantennas = 10; 
  combine_npad = 3; 
  combine_unfiltered = true; 

  saturation_threshold = 2500; 

  start_pol = AnitaPol::kHorizontal; 
  end_pol = AnitaPol::kVertical; 

  peak_isolation_requirement = 20; 

  fine_peak_finding_option = FinePeakFindingQuadraticFit25; 
  nmaxima = 2; 
  use_bin_center = false; 
  scale_by_cos_theta = false; 

  bw_ndb = 6; 
  spectral_fit_start = 0.22; 
  spectral_fit_stop =0.7; 

  noise_estimate_t0 = 70; 
  noise_estimate_t1 = 100; 

  response_option = ResponseNone; 
  response_npad = 50; 


  if (config) loadFromFile(config);

  deconvolution_method = &kDefaultDeconvolution; 
}


const char * UCorrelator::AnalysisConfig::getPeakFindingString(FinePeakFindingOption_t opt) 
{
  return peakfinders[opt]; 
}

const char * UCorrelator::AnalysisConfig::getResponseString(ResponseOption_t opt) 
{
  return responses[opt]; 
}



