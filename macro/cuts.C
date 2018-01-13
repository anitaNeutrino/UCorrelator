#include "FFTtools.h" 
#define CUTS_LOADED 
//neutrino search, every cuts is specialized for V pol .
TCut isReal ( "peak[1][0].value > 0 && peak[1][1].value > 0 && flags.isRF==1 && peak[1][0].theta < 60 &&  peak[1][0].theta > -50 && flags.hasGlitch == 0"); 
TCut isWaisH( "flags.pulser == 1 && flags.isRF == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[0][0].theta - wais.theta,360,0)) < 2 && flags.hasGlitch == 0"); 
TCut isWaisV( "flags.pulser == 5 && flags.isRF == 1 && abs(FFTtools::wrap(peak[1][0].phi - wais.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[1][0].theta - wais.theta,360,0)) < 2 && flags.hasGlitch == 0"); 
TCut notWais( "flags.pulser != 1 && flags.pulser != 5"); 
TCut isSun("flags.pulser == 0 && flags.isRF == 1 && abs(FFTtools::wrap(peak[1][0].phi - sun.phi,360,0)) < 20"); 
TCut isNorth("abs(FFTtools::wrap(peak[1][0].phi - heading,360,0)) < 90"); 
TCut triggered(" abs(peak[1][0].hwAngle < 60)"); 
TCut notMasked ("!peak[1][0].masked || !peak[1][0].masked_xpol"); 
TCut aboveHorizon (" peak[1][0].theta < 0"); 
// TCut isntSidelobe (" peak[1][0].phi_separation > 10"); 
TCut isGlitch("flags.hasGlitch == 1");
TCut badReconstructionCut = " abs(peak[1][0].dphi_rough) < 4 && abs(peak[1][0].dtheta_rough) <3"; 
// TCut notBlast = "flags.meanPower[3]/flags.meanPower[1]>0.2 && flags.meanPower[3]/flags.meanPower[1]<1.7 && flags.maxBottomToTopRatio[0] >0.9 && flags.maxBottomToTopRatio[0] < 2.6 && flags.maxBottomToTopRatio[1] >0.9 && flags.maxBottomToTopRatio[1] < 2.6";
TCut notBlast = "flags.isPayloadBlast==0";
TCut notTooFiltered = "flags.meanPowerFiltered[0] / flags.meanPower[0] > 0.2"; 
TCut thermal_sample = isReal && notWais && badReconstructionCut && notBlast && triggered && notMasked; 
TCut isMC = "abs(FFTtools::wrap(peak[1][0].phi-mc.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[1][0].theta-mc.theta,360,0)) < 3 && ( ( mc.wf[0].peakHilbert > mc.wf[1].peakHilbert  &&  Iteration$ < 5) || ( mc.wf[1].peakHilbert > mc.wf[0].peakHilbert && Iteration$ >= 5))"; 
TCut anyMC = "Sum$((abs(FFTtools::wrap(peak[1][0].phi-mc.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[1][0].theta-mc.theta,360,0)) < 3)  > 0)"; 