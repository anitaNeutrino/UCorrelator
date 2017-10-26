#ifndef _UCORRELATOR_DELTA_T_H
#define _UCORRELATOR_DELTA_T_H

#include "AntennaPositions.h" 
#include "AnitaGeomTool.h"
#include "FFTtools.h"
#include "TrigCache.h" 
#include <TMVA/Reader.h>


#ifndef DEG2RAD
#define DEG2RAD (M_PI/180.); 
#endif

namespace UCorrelator
{

  /** Placeholder for now ...  this should be updated for A3*/ 
  inline double getAntennaGroupDelay(double phidiff, double theta) 
  {
    theta-=10;
    const double c1 = 1.45676e-8; 
    const double c2 = 5.01452e-6; 

    double totalAngle2=theta*theta+phidiff*phidiff;//TODO
    if (totalAngle2>50*50) return (50*50*50*50 * c1 - 50*50 * c2); 

    double delayTime=totalAngle2 * totalAngle2 * c1 - totalAngle2 * c2; 

    return delayTime;
  }


  /**Geometric delay between the phase center of an antenna and a center point  */
  inline double getDeltaTtoCenter(int ant1, double phi, double theta, AnitaPol::AnitaPol_t pol, bool includeGroupDelay = false) 
  {
    double th = theta * DEG2RAD; 
    const AntennaPositions * ap = AntennaPositions::instance(); 
    double ph1_deg = (phi- ap->phiAnt[pol][ant1]) ; 
    double ph1  = ph1_deg * DEG2RAD; 

    double part1=ap->zAnt[pol][ant1]*tan(th) - ap->rAnt[pol][ant1] * cos(ph1);
    
    double geomDelay=1e9*((cos(th) * part1)/C_LIGHT);    //returns time in ns
    

    if (includeGroupDelay)
    {
      geomDelay +=  getAntennaGroupDelay(FFTtools::wrap(ph1_deg,360,0), theta);
    }

    return geomDelay;
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

  /** Geometric delay between antennas using cached trig values 
   *
   *  TODO: This has to be vectorized somehow, at least over one dimension, as it's currently one of the bottlenecks. 
   *
   * */ 
  inline double getDeltaTFast(int ant1, int ant2, int phibin, int thetabin, AnitaPol::AnitaPol_t pol, const TrigCache * cache, bool includeGroupDelay = false) 
  {
    const AntennaPositions * ap = cache->ap; 
    const int nphi = cache->nphi; 
    double part1=ap->zAnt[pol][ant1]*cache->tan_theta[thetabin] - ap->rAnt[pol][ant1] * cache->cos_phi[2*(phibin + nphi * ant1) + pol];
    double part2=ap->zAnt[pol][ant2]*cache->tan_theta[thetabin] - ap->rAnt[pol][ant2] * cache->cos_phi[2*(phibin + nphi * ant2) + pol];
    
    double geomDelay=(1.e9/C_LIGHT)*(cache->cos_theta[thetabin] * (part1 - part2));    //returns time in ns


    if (includeGroupDelay)
    {
      double ph1_deg = FFTtools::wrap((cache->phi[phibin]- ap->phiAnt[pol][ant1]),360,0) ; 
      double ph2_deg = FFTtools::wrap((cache->phi[phibin]- ap->phiAnt[pol][ant2]),360,0) ; 
      double theta_deg = cache->theta[thetabin];
      geomDelay +=  getAntennaGroupDelay(ph1_deg, theta_deg) - getAntennaGroupDelay(ph2_deg, theta_deg); 
    }

    return geomDelay;
  }


    //get the deltaT given (ant1, ant2, pol, theta ,phi) from Artificial Neurual Network(ANN).
  inline double getDeltaTFromANN(int ant1, int ant2, int phibin, int thetabin, AnitaPol::AnitaPol_t pol,const UCorrelator::TrigCache * cache , TMVA::Reader* reader){
    
    if(ant1 > ant2){
      int tmp = ant1;
      ant1 = ant2;
      ant2 = tmp;
    }
    char antPair[100];
    sprintf (antPair, "antPair_%d_%d_%d", ant1, ant2, pol);
    float phi_rad, theta_rad;

    // two input variable theta and phi to the ANN.
    reader->AddVariable("phi",&phi_rad);
    reader->AddVariable("theta",&theta_rad);

    phi_rad = cache->phi[phibin]*DEG2RAD;
    theta_rad = cache->theta[thetabin]*DEG2RAD;

    //get the results from ANN and return
    return reader->EvaluateRegression(antPair)[0];
  }




  
}

#endif
