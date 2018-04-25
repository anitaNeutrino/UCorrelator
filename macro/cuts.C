#include "FFTtools.h" 
#define CUTS_LOADED 
//neutrino search, every cuts is specialized for V pol .
TCut isReal ( "flags.isRF==1 && mostImpulsivePeak(2).theta < 60 &&  mostImpulsivePeak(2).theta > -50"); 
// TCut allWais("flags.pulser == 1 || flags.pulser == 5");
TCut allWais("flags.pulser == 5 || flags.pulser == 1");
TCut isWaisH( "flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[0][0].theta - wais.theta,360,0)) < 2"); 
TCut isWaisV( "flags.pulser == 5 && abs(FFTtools::wrap(peak[1][0].phi - wais.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[1][0].theta - wais.theta,360,0)) < 2"); 
TCut notWais( "flags.pulser != 1 && flags.pulser != 5"); 
TCut isSun("flags.pulser == 0 && flags.isRF == 1 && abs(FFTtools::wrap(peak[1][0].phi - sun.phi,360,0)) < 20"); 
TCut isNorth("abs(FFTtools::wrap(peak[1][0].phi - heading,360,0)) < 90");
TCut notBadReconstruction(" abs(mostImpulsivePeak(2).dphi_rough) < 6 && abs(mostImpulsivePeak(2).dtheta_rough) <3");  
TCut triggered(" abs(mostImpulsivePeak(2).hwAngle) < 60"); 
TCut notMasked ("!mostImpulsivePeak(2).masked"); 
TCut aboveHorizontal("mostImpulsivePeak(2).theta<0");
TCut belowHorizon("mostImpulsivePeak(2).theta>6");
TCut betweenHorizonAndHorizontal("mostImpulsivePeak(2).theta>0 && mostImpulsivePeak(2).theta<5");
// TCut isntSidelobe (" peak[1][0].phi_separation > 10");
TCut notGlitch("flags.hasGlitch == 0"); 
TCut isGlitch("flags.hasGlitch == 1");
// TCut notBlast = "flags.meanPower[3]/flags.meanPower[1]>0.2 && flags.meanPower[3]/flags.meanPower[1]<1.7 && flags.maxBottomToTopRatio[0] >0.9 && flags.maxBottomToTopRatio[0] < 2.6 && flags.maxBottomToTopRatio[1] >0.9 && flags.maxBottomToTopRatio[1] < 2.6";
TCut notBlast("flags.isPayloadBlast==0");
TCut notTooFiltered("flags.meanPowerFiltered[0] / flags.meanPower[0] > 0.2"); 
TCut notStrongCW("mostImpulsiveCoherent(2).snr/mostImpulsiveDeconvolved(2).snr < 5.5");
TCut snrCut("mostImpulsiveDeconvolvedFiltered(2).snr>9.5");
// TCut isMC = "abs(FFTtools::wrap(peak[1][0].phi-mc.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[1][0].theta-mc.theta,360,0)) < 2 && ( ( mc.wf[0].peakHilbert > mc.wf[1].peakHilbert  &&  Iteration$ < 5) || ( mc.wf[1].peakHilbert > mc.wf[0].peakHilbert && Iteration$ >= 5))"; 
// TCut isMC = "mc.wf[1].peakHilbert > mc.wf[0].peakHilbert"; 
TCut goodPointingMC("abs(FFTtools::wrap(peak[1][0].phi-mc.phi,360,0) < 5 && abs(FFTtools::wrap(peak[1][0].theta-mc.theta,360,0)) < 2)");
TCut notHical0("Hical2::isHical(eventNumber, header->triggerTime, FFTtools::wrap(anitaLocation.heading - peak[0][0].phi, 360))!=1");
TCut notHical1("Hical2::isHical(eventNumber, header->triggerTime, FFTtools::wrap(anitaLocation.heading - peak[0][1].phi, 360))!=1");
TCut notHical2("Hical2::isHical(eventNumber, header->triggerTime, FFTtools::wrap(anitaLocation.heading - peak[0][2].phi, 360))!=1");
TCut notHical = notHical0 && notHical1 && notHical2;
TCut qualityCut = isReal && notGlitch && notBadReconstruction && notBlast && triggered && notMasked && notStrongCW && notHical;
TCut goodPointingWais = (isWaisH || isWaisV);
TCut thermal_sample = qualityCut;
TCut wais_sample =  qualityCut && goodPointingWais;
TCut mc_sample = qualityCut && goodPointingMC;
