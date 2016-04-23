#include "AnalysisConfig.h" 


const int wais_hpol_time_offset = 93; 
const int wais_vpol_time_offset = -99757; 
const int ldb_hpol_time_offset = 0;  //TODO
const int ldb_vpol_time_offset = 0; 

UCorrelator::AnalysisConfig::AnalysisConfig() 
  : 
    wais_hpol(wais_hpol_time_offset, 800e3, 1e3), 
    wais_vpol(wais_vpol_time_offset, 800e3, 1e3), 
    ldb_hpol(ldb_hpol_time_offset, 800e3, 1e3), 
    ldb_vpol(ldb_vpol_time_offset, 800e3, 1e3) 
{
  correlator_nphi = 180; 
  correlator_ntheta = 90; 
  correlator_theta_lowest = 60; 
  correlator_theta_highest = -25; 


  zoomed_nphi = 50; 
  zoomed_ntheta = 50; 
  zoomed_dphi = 0.3; 
  zoomed_dtheta = 0.3; 


  combine_nantennas = 10; 
  combine_npad = 3; 

  saturation_threshold = 2500; 

  start_pol = AnitaPol::kHorizontal; 
  end_pol = AnitaPol::kVertical; 

  peak_isolation_requirement = 10; 

  fine_peak_finding_option = FinePeakFindingQuadraticFit16; 
  nmaxima = 2; 
  use_bin_center = false; 

  bw_ndb = 3; 

  noise_estimate_t0 = 70; 
  noise_estimate_t1 = 100; 

}

