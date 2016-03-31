#include "AntennaPositions.h" 
#include "AnitaGeomTool.h"
#include "assert.h"
#include "FFTtools.h"
#include <map>


UCorrelator::AntennaPositions * UCorrelator::AntennaPositions::theInstance = 0; 

#ifndef RAD2DEG
#define RAD2DEG (180/M_PI)
#endif




UCorrelator::AntennaPositions::AntennaPositions()
{

  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {
    for (int ant  = 0; ant < NUM_SEAVEYS; ant++) 
    {
//      printf("%d %d\n",pol,ant); 
      AnitaGeomTool * geom = AnitaGeomTool::Instance(); 
      phiAnt[pol][ant] = geom->getAntPhiPositionRelToAftFore(ant,(AnitaPol::AnitaPol_t)pol) * RAD2DEG; 
      rAnt[pol][ant] = geom->getAntR(ant); 
      zAnt[pol][ant] = geom->getAntZ(ant); 
    }
  }
}



void UCorrelator::AntennaPositions::getClosestAntennas(double phi, int N, int * closest, uint64_t disallowed ) const
{

  assert(N < NUM_SEAVEYS); 

  std::multimap<double,int> dphis; 
  for (int i = 0; i < NUM_SEAVEYS; i++) 
  {

    double dphi = disallowed & (1 << i) ? 360 : fabs(FFTtools::wrap(phiAnt[0][i], 360, 0)); 
    dphis.insert(std::pair<double,int>(dphi,i)); 
  }


  int Nused = 0; 

  for (std::multimap<double,int>::const_iterator it = dphis.begin(); it !=dphis.end(); it++) 
  {
    closest[Nused++] = (*it).second; 
    if (Nused == N) break; 
  }
}



