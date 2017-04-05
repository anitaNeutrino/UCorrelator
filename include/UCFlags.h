#ifndef _UCORRELATOR_FLAGS_H
#define _UCORRELATOR_FLAGS_H

/** \file This file is full of a bunch of methods that check various things about events */ 

#include <stdint.h>
class UsefulAnitaEvent; 

namespace UCorrelator
{
  namespace flags 
  {
    /** Checks for missing antennas and marks the bitmasks */ 
    int checkEmpty(const UsefulAnitaEvent *ev,  uint64_t *hempty, uint64_t *vempty); 


  }
}


#endif
