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

  /*=====
    TH1::GetRMS(3) stopped working so here is a function to just calculate it by hand
   */
  double getZRMS(const TH2*);



  namespace image
  {

    TH1 * getPctileProjection(const TH2 * h, int axis = 1, double pct = 0.5, bool ignoreEmpty = true); 

  }; 
}; 

#endif
