#include "UCFlags.h" 
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h" 




int UCorrelator::flags::checkEmpty(const UsefulAnitaEvent *ev, ULong64_t * save_h, ULong64_t * save_v)
{
  

  int nmissing = 0; 
  for (int i = 0; i < NUM_SEAVEYS; i++) 
  {
      int hindex = AnitaGeomTool::getChanIndexFromAntPol(i,AnitaPol::kHorizontal); 
      int vindex = AnitaGeomTool::getChanIndexFromAntPol(i,AnitaPol::kVertical); 
      const double *yh = ev->fVolts[hindex]; 
      const double *yv = ev->fVolts[vindex]; 

      bool hbad = true; 
      for (int j = 0; j < ev->fNumPoints[hindex]; j++)
      {
        if (yh[j]) 
        {
          hbad = false; 
          break; 
        }
      }

      bool vbad = true; 
      for (int j = 0; j < ev->fNumPoints[vindex]; j++)
      {
        if (yv[j])
        {
          vbad = false; 
          break; 
        }
      }


      if (hbad && save_h) *save_h |= (1 << i); 
      if (vbad && save_v) *save_v |= (1 << i); 


      if (hbad) nmissing++; 
      if (vbad) nmissing++; 

  }


  return nmissing; 
}
