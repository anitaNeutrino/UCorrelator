#include "FFTtools.h" 
#define CUTS_LOADED 

TCut isReal ( "peak[][].value > 0 && flags.isRF==1 && peak[][].theta < 60 && peak[][].theta > -50"); 
TCut isWais( "flags.pulser == 1 && flags.isRF == 1 && abs(FFTtools::wrap(peak[][].phi - wais.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[][].theta - wais.theta,360,0)) < 5"); 
TCut isSun("flags.pulser == 0 && flags.isRF == 1 && abs(FFTtools::wrap(peak[][].phi - sun.phi,360,0)) < 20"); 
TCut isNorth("abs(FFTtools::wrap(peak[][].phi - heading,360,0)) < 90"); 
TCut triggered(" abs(peak[][].hwAngle < 60)"); 
TCut notMasked ("!peak[][].masked || !peak[][].masked_xpol"); 
TCut aboveHorizon (" peak[][].theta < 0"); 
TCut isntSidelobe (" peak[][].phi_separation > 10"); 
TCut antiglitch ("coherent[][].peakHilbert / coherent[][].peakVal > 1 && coherent[][].peakHilbert / coherent[][].peakVal < 2.5"); 
TCut badReconstructionCut = " abs(peak.dphi_rough) < 4 && abs(peak.dtheta_rough) <3"; 
TCut blastCut = "!flags.isPayloadBlast && flags.maxBottomToTopRatio < 2.6 && flags.maxBottomToTopRatio > 1";
TCut brightestPeak = "(Iteration$ == 0 && peak[0][0].value > peak[1][0].value) || Iteration$ == 5 && (peak[1][0].value > peak[0][0].value)";
TCut notTooFiltered = "flags.meanPowerFiltered[0] / flags.meanPower[0] > 0.2"; 
TCut not460 = "abs(coherent[][].peakFrequency[0]-0.46) > 0.02"; 
TCut thermal_sample = isReal && brightestPeak && !isWais && !isSun && !isNorth && aboveHorizon && isntSidelobe && badReconstructionCut && blastCut && triggered; 

