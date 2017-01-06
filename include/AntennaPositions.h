#ifndef _UCORRELATOR_ANTENNA_POSITIONS_H
#define _UCORRELATOR_ANTENNA_POSITIONS_H

#include "AnitaConventions.h" 
#include "AnitaVersion.h" 
#include "stdint.h"

namespace UCorrelator
{
  /** This class keeps the positions of the antennas, and has some related methods. It is a singleton. */ 
  class AntennaPositions
  {

    static const AntennaPositions * instances[NUM_ANITAS+1]; 

    AntennaPositions(int v); 

    public: 
      /** Retrieve an instance */ 
      static const AntennaPositions * instance(int version = 0)
      { 

        if (__builtin_expect(!version, 1)) 
        {
          version = AnitaVersion::get();  
        }

        if (__builtin_expect(!instances[version],0)) 
        {
          instances[version] = new AntennaPositions(version); 
        }

        return instances[version];
      }

      /** Find closest N antennas to phi. Results put into closest, which should have sufficient room. Disallowed is a bitmap of antenna numbers that should be excluded. Returns number found (could be less than number requested if too many disallowed)*/
      int getClosestAntennas(double phi, int N, int * closest, uint64_t disallowed = 0) const; 

      /** antenna phi positions (degrees)*/
      double phiAnt[2][NUM_SEAVEYS]; 

      /** antenna r positions (m) */
      double rAnt[2][NUM_SEAVEYS]; 

      /** antenna z positions (m) */
      double zAnt[2][NUM_SEAVEYS]; 
  }; 
}


#endif
