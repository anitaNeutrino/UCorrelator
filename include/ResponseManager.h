#ifndef UCORRELATOR_RESPONSE_MANAGER_H 
#define UCORRELATOR_RESPONSE_MANAGER_H

/* This class does the dirty work of loading appopriate responses
 * for each antenna. See README.response for more details.
 *
 *
 *   Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 * */ 

#include "AnalysisConfig.h" 
#include <vector> 

namespace UCorrelator
{
  class AbstractResponse; 
  class ResponseManager
  {

    public: 

      ResponseManager(AnalysisConfig::ResponseOption_t option); 
      ResponseManager(const char * responseDir); 

      const AbstractResponse * response(int pol, int iant) const { return responses[iant][pol]; } 

      virtual ~ResponseManager(); 

      
     private: 
     int loadResponsesFromDir(const char * dir); 
     const AbstractResponse* responses[NUM_SEAVEYS][2]; 
     std::vector<AbstractResponse*> response_store; 


  };


}


#endif

