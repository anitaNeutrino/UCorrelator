#ifndef UCORRELATOR_SHAPE_PARAMETERS_HH
#define UCORRELATOR_SHAPE_PARAMETERS_HH

class TGraph; 

namespace UCorrelator
{
  namespace shape
  {

    double getRiseTime(const TGraph * g, double val_min, double val_max); 
    double getFallTime(const TGraph * g, double val_min, double val_max); 
    double getWidth(const TGraph * g, double val_thresh); 
  }

}



#endif 
