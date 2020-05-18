#ifndef _UC_CACHED_FC_H
#define _UC_CACHED_FC_H

#include "TFeldmanCousins.h" 


namespace UCorrelator

{
class CachedFC
{

  public: 
    CachedFC(double cl, int max_sig, double max_bg, double dbg = 0.1);

    virtual ~CachedFC(); 
    double upperLimit(int nsig, double bg); 

  private: 
    TFeldmanCousins fc; 

    //hide this 
    void* interp; 
}; 

}

#endif
