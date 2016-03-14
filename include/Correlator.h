#ifndef UCORRELATOR_CORRELATOR_H
#define UCORRELATOR_CORRELATOR_H


class FilteredAnitaEvent; 
class TH2; 

namespace UCorrelator
{
  class Correlator
  {
    public:
      Correlator(int nphi, double phimin, double phimax, int ntheta, double theta_min, double theta_max); 


      void compute(FilteredAnitaEvent * event); 

      TH2 * getHist(); 


  }; 



}



#endif
