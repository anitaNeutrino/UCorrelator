#ifndef _UCORRELATOR_FLAGS_H
#define _UCORRELATOR_FLAGS_H

/** This file is full of a bunch of methods that check various things about events */ 

#include <stdint.h>
class UsefulAnitaEvent; 

namespace UCorrelator
{
  namespace flags 
  {
    /** Checks event for saturation (any |value| > threshold ), returning the total number of saturated channels
     * Returns number of saturated traces and optionally populates bits of hsta and vsat (hopefully Anita-9000 doesn't have more than 64 antennas!) 
     */
    int checkSaturation(const UsefulAnitaEvent *ev, uint64_t *hsat = 0, uint64_t * vsat = 0, double threshold_in_mV = 1500); 



  }
}


#endif
