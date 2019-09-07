#ifndef UCORRELATOR_SHAPE_PARAMETERS_HH
#define UCORRELATOR_SHAPE_PARAMETERS_HH

class TGraph; 

namespace UCorrelator
{
  namespace shape
  {

    double getRiseTime(const TGraph * g, double val_min, double val_max, int peakBin = -1); 
    double getFallTime(const TGraph * g, double val_min, double val_max, int peakBin = -1); 
    /* Gets the width around the maximum. if PeakBin is >=0, that is assumed to be the peak (otherwise it's calculated)*/  
    double getWidth(const TGraph * g, double val_thresh, int * istart = 0, int * iend = 0, int peakBin = -1); 
  }

}



#endif 
