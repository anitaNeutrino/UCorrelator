#ifndef _UCORRELATOR_DELTA_T_H
#define _UCORRELATOR_DELTA_T_H

#include "AntennaPositions.h" 
#include "AnitaGeomTool.h"
#include "FFTtools.h"
#include "TrigCache.h" 


#ifndef DEG2RAD
#define DEG2RAD (M_PI/180.); 
#endif

namespace UCorrelator
{
  /** Placeholder for now ...  this should be updated for A3*/ 
  inline double getAntennaGroupDelay(double phidiff, double theta) 
  {
    theta-=10;
    double totalAngle=sqrt(theta*theta+phidiff*phidiff);//TODO
    if (totalAngle>50) totalAngle=50;

    double delayTime=(totalAngle*totalAngle*totalAngle*totalAngle)*1.45676e-8;
    delayTime-=(totalAngle*totalAngle)*5.01452e-6;

    return delayTime;
  }


  /**Geometric delay between antennas  */
  inline double getDeltaT(int ant1, int ant2, double phi, double theta, AnitaPol::AnitaPol_t pol, bool includeGroupDelay = false) 
  {
    double th = theta * DEG2RAD; 
    const AntennaPositions * ap = AntennaPositions::instance(); 
    double ph1_deg = (phi- ap->phiAnt[pol][ant1]) ; 
    double ph2_deg = (phi- ap->phiAnt[pol][ant2]) ; 
    double ph1  = ph1_deg * DEG2RAD; 
    double ph2  = ph2_deg * DEG2RAD; 

    double part1=ap->zAnt[pol][ant1]*tan(th) - ap->rAnt[pol][ant1] * cos(ph1);
    double part2=ap->zAnt[pol][ant2]*tan(th) - ap->rAnt[pol][ant2] * cos(ph2); 
    
    double geomDelay=1e9*((cos(th) * (part1 - part2))/C_LIGHT);    //returns time in ns
    

    if (includeGroupDelay)
    {
      geomDelay +=  getAntennaGroupDelay(FFTtools::wrap(ph1_deg,360,0), theta) - getAntennaGroupDelay(FFTtools::wrap(ph2_deg,360,0), theta); 
    }

    return geomDelay;
  }

  /** Geometric delay between antennas using cached trig values */ 
  inline double getDeltaTFast(int ant1, int ant2, int phibin, int thetabin, AnitaPol::AnitaPol_t pol, const TrigCache * cache, bool includeGroupDelay = false) 
  {
    const AntennaPositions * ap = AntennaPositions::instance(); 
    const int nphi = cache->nphi; 
    double part1=ap->zAnt[pol][ant1]*cache->tan_theta[thetabin] - ap->rAnt[pol][ant1] * cache->cos_phi[2*(phibin + nphi * ant1) + pol];
    double part2=ap->zAnt[pol][ant2]*cache->tan_theta[thetabin] - ap->rAnt[pol][ant2] * cache->cos_phi[2*(phibin + nphi * ant2) + pol];
    
    double geomDelay=1e9*((cache->cos_theta[thetabin] * (part1 - part2))/C_LIGHT);    //returns time in ns


    if (includeGroupDelay)
    {
      double ph1_deg = FFTtools::wrap((cache->phi[phibin]- ap->phiAnt[pol][ant1]),360,0) ; 
      double ph2_deg = FFTtools::wrap((cache->phi[phibin]- ap->phiAnt[pol][ant2]),360,0) ; 
      double theta_deg = cache->theta[thetabin];
      geomDelay +=  getAntennaGroupDelay(ph1_deg, theta_deg) - getAntennaGroupDelay(ph2_deg, theta_deg); 
    }

    return geomDelay;
  }

  
}

#endif
