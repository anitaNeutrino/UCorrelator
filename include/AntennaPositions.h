#ifndef _UCORRELATOR_ANTENNA_POSITIONS_H
#define _UCORRELATOR_ANTENNA_POSITIONS_H

#include "AnitaConventions.h" 

namespace UCorrelator
{
  class AntennaPositions
  {

    static AntennaPositions * theInstance; 
    public: 
      static AntennaPositions * instance(); 

      double phiAntByTrace[NUM_DIGITZED_CHANNELS]; 
      double rAntByTrace[NUM_DIGITZED_CHANNELS]; 
      double zAntByTrace[NUM_DIGITZED_CHANNELS]; 

      int antByTrace[NUM_DIGITZED_CHANNELS]; 
      AnitaPol::AnitaPol_t polByTrace[NUM_DIGITZED_CHANNELS]; 
  }; 
}


#endif
