#include "AntennaPositions.h" 
#include "AnitaGeomTool.h"
#include "assert.h"
#include "FFTtools.h"
#include <map>



#ifndef RAD2DEG
#define RAD2DEG (180/M_PI)
#endif



const UCorrelator::AntennaPositions * UCorrelator::AntennaPositions::instances[NUM_ANITAS+1] = {0}; 


UCorrelator::AntennaPositions::AntennaPositions(int version)
{

  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {
    for (int ant  = 0; ant < NUM_SEAVEYS; ant++) 
    {
//      printf("%d %d\n",pol,ant); 
      AnitaGeomTool * geom = AnitaGeomTool::Instance(version); 
      phiAnt[pol][ant] = geom->getAntPhiPositionRelToAftFore(ant,(AnitaPol::AnitaPol_t)pol) * RAD2DEG; 
      rAnt[pol][ant] = geom->getAntR(ant,(AnitaPol::AnitaPol_t) pol); 
      zAnt[pol][ant] = geom->getAntZ(ant,(AnitaPol::AnitaPol_t) pol); 
    }
  }

  v = version; 
}



int UCorrelator::AntennaPositions::getClosestAntennas(double phi, int N, int * closest, uint64_t disallowed , AnitaPol::AnitaPol_t pol) const
{

  assert(N < NUM_SEAVEYS); 

  int pol_ind = pol == AnitaPol::kHorizontal ? 0 : 1; 
  std::multimap<double,int> dphis; 
  for (int i = 0; i < NUM_SEAVEYS; i++) 
  {

    double dphi = disallowed & (1 << i) ? 360 : fabs(FFTtools::wrap(phi-phiAnt[pol_ind][i], 360, 0)); 
    dphis.insert(std::pair<double,int>(dphi,i)); 
  }

  int Nused = 0; 

//  printf("Closest antennas to %f: ",phi); 
  for (std::multimap<double,int>::const_iterator it = dphis.begin(); it !=dphis.end(); it++) 
  {
//    printf("  %d %f\n", (*it).second, (*it).first); 
    
    if ((*it).second >=360) break; 
    closest[Nused++] = (*it).second; 
    if (Nused == N) break; 
  }
//  printf("\n"); 
  return Nused; 
}



double UCorrelator::AntennaPositions::distance(int i, int j, AnitaPol::AnitaPol_t pol) const
{

  AnitaGeomTool * geom = AnitaGeomTool::Instance(v); 

  double x0,y0,z0; 
  double x1,y1,z1; 
  geom->getAntXYZ(i,x0,y0,z0, pol); 
  geom->getAntXYZ(j,x1,y1,z1, pol); 


  double dx = x0-x1; 
  double dy = y0-y1; 
  double dz = z0-z1; 
  return sqrt( dx*dx + dy*dy + dz*dz); 

}

