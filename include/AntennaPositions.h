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
      static AntennaPositions * instance(); 

      void getClosestAntennas(double phi, int N, int * closest, uint64_t disallowed = 0); 
      double phiAntByTrace[NUM_DIGITZED_CHANNELS]; 
      double phiAntByPolAnt[2][NUM_DIGITZED_CHANNELS/2]; 
      double rAntByTrace[NUM_DIGITZED_CHANNELS]; 
      double rAntByPolAnt[2][NUM_DIGITZED_CHANNELS/2]; 
      double zAntByTrace[NUM_DIGITZED_CHANNELS]; 
      double zAntByPolAnt[2][NUM_DIGITZED_CHANNELS/2]; 
      int traceByPolAnt[2][NUM_DIGITZED_CHANNELS/2]; 
      int antByTrace[NUM_DIGITZED_CHANNELS]; 
      AnitaPol::AnitaPol_t polByTrace[NUM_DIGITZED_CHANNELS]; 
  }; 
}


#endif
