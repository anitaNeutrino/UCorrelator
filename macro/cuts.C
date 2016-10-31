{

#include "FFTtools.h" 

  TCut isReal ( "peak[][].value > 0 && flags.isRF==1 && peak[][].theta < 60 && peak[][].theta > -60"); 
  TCut isWais( "flags.pulser == 1 && flags.isRF == 1 && abs(FFTtools::wrap(peak[][].phi - wais.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[][].theta - wais.theta,360,0)) < 5"); 
  TCut isSun("flags.pulser == 0 && flags.isRF == 1 && abs(FFTtools::wrap(peak[][].phi - sun.phi,360,0)) < 20"); 
  TCut isNorth("abs(FFTtools::wrap(peak[][].phi - heading,360,0)) < 90"); 
  TCut triggered(" abs(peak[][].triggerAngle < 60"); 
  TCut aboveHorizon (" peak[][].theta < 0"); 
  TCut isntSidelobe (" peak[][].phi_separation > 5"); 
  TCut badReconstructionCut = " abs(peak.dphi_rough) < 5 && abs(peak.dtheta_rough) <5"; 
  TCut blastCut = "flags.meanPowerFiltered[0] < 1e6 && flags.maxBottomToTopRatio[0] && flags.maxBottomToTopRatio[1] < 5 "; 

  TCut thermal_sample = isReal && !isWais && !isSun && !isNorth && aboveHorizon && isntSidelobe && badReconstructionCut && blastCut; 

}
