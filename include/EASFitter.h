#ifndef UCORRELATOR_EAS_FITTER_H
#define UCORRELATOR_EAS_FITTER_H

#include <vector> 

class TCanvas; 
class FilteredAnitaEvent; 
class TObject; 

namespace AnitaResponse
{
  class ResponseManager; 
}
namespace UCorrelator
{


  struct EASFitResult
  {
    struct MinimizationResult 
    {
      int status; //minimization status 
      double chisq; 
      double initialA;
      double A; 
      double sigmaA; 
      double initialT; 
      double T; 
      double sigmaT; 
    }; 

    std::vector<MinimizationResult> result; 
  };


  class EASFitter 
  {
    public:

      EASFitter(const AnitaResponse::ResponseManager * responses); 

      ~EASFitter(); 


      int fitEvent(int nants, const int * ants, 
                   const FilteredAnitaEvent * event,
                   double phi, double theta, bool dedisperse); 

      unsigned nResults() const { return results.size(); } 
      const EASFitResult * result(int i) const { return &results.at(i); } 
      TCanvas * plot(int i) { return plots.at(i); } 

    private: 
      const AnitaResponse::ResponseManager * rm; 
      std::vector<EASFitResult> results; 
      std::vector<TCanvas *> plots; 
      std::vector<TObject*> save; 
  }; 



} 

#endif
