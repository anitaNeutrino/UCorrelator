#ifndef _UCORRELATOR_DELTA_T_H
#define _UCORRELATOR_DELTA_T_H

#include "AntennaPositions.h" 
#include "AnitaGeomTool.h"
#include "TrigCache.h" 


#ifndef DEG2RAD
#define DEG2RAD (M_PI/180); 
#endif

namespace UCorrelator
{

  //Geometric delay between antennas; 
  inline double getDeltaT(int ant1, int ant2, double phi, double theta, AnitaPol::AnitaPol_t pol) 
  {
    double th = theta * DEG2RAD; 
    const AntennaPositions * ap = AntennaPositions::instance(); 
    double ph1 = DEG2RAD * (phi  -ap->phiAnt[pol][ant1]) ; 
    double ph2 = DEG2RAD * (phi - ap->phiAnt[pol][ant2]) ; 
    double part1=ap->zAnt[pol][ant1]*tan(th) - ap->rAnt[pol][ant1] * cos(ph1);
    double part2=ap->zAnt[pol][ant2]*tan(th) - ap->rAnt[pol][ant2] * cos(ph2); 
    
    double geomDelay=1e9*((cos(th) * (part1 - part2))/C_LIGHT);    //returns time in ns
    
    return geomDelay;
  }

  //geometric delay using cache
  inline double getDeltaTFast(int ant1, int ant2, int phibin, int thetabin, AnitaPol::AnitaPol_t pol, const TrigCache * cache) 
  {
    const AntennaPositions * ap = AntennaPositions::instance(); 
    const int nphi = cache->nphi; 
    double part1=ap->zAnt[pol][ant1]*cache->tan_theta[thetabin] - ap->rAnt[pol][ant1] * cache->cos_phi[2*(phibin + nphi * ant1) + pol];
    double part2=ap->zAnt[pol][ant2]*cache->tan_theta[thetabin] - ap->rAnt[pol][ant2] * cache->cos_phi[2*(phibin + nphi * ant2) + pol];
    
    double geomDelay=1e9*((cache->cos_theta[thetabin] * (part1 - part2))/C_LIGHT);    //returns time in ns

    return geomDelay;
  }

  
  /** Placeholder for now ... */ 
  inline double getAntennaGroupDelay(double phidiff, double theta) 
  {
    theta-=10;
    double totalAngle=sqrt(theta*theta+phidiff*phidiff);//TODO
    if (totalAngle>50) totalAngle=50;

    double delayTime=(totalAngle*totalAngle*totalAngle*totalAngle)*1.45676e-8;
    delayTime-=(totalAngle*totalAngle)*5.01452e-6;

    return delayTime;
  }

}

#endif
