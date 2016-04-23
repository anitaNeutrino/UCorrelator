#ifndef _UCORRELATOR_ANTENNA_POSITIONS_H
#define _UCORRELATOR_ANTENNA_POSITIONS_H

#include "AnitaConventions.h" 
#include "stdint.h"

namespace UCorrelator
{
  class AntennaPositions
  {

    static AntennaPositions * theInstance; 
    AntennaPositions(); 

    public: 
      static const AntennaPositions * instance() { if (!theInstance) theInstance = new AntennaPositions ; return theInstance; }

      void getClosestAntennas(double phi, int N, int * closest, uint64_t disallowed = 0) const; 

      double phiAnt[2][NUM_SEAVEYS]; 
      double rAnt[2][NUM_SEAVEYS]; 
      double zAnt[2][NUM_SEAVEYS]; 
  }; 
}


#endif
