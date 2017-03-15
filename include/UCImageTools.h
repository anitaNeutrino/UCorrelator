#ifndef UCORRELATOR_IMAGE_TOOLS_H
#define UCORRELATOR_IMAGE_TOOLS_H

/** Image tools for working with TH2*'s
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 **/


class TH1; 
class TH2; 


namespace UCorrelator
{


  namespace image
  {

    TH1 * getPctileProjection(const TH2 * h, int axis = 1, double pct = 0.5, bool ignoreEmpty = true); 

  }; 
}; 

#endif
