#ifndef ANITA_NOISE_MACHINE
#define ANITA_NOISE_MACHINE

#include "AnitaConventions.h"
#include "AnitaNoiseSummary.h"
#include "AnitaEventSummary.h"
#include "FilteredAnitaEvent.h"
#include "Analyzer.h"
#include "UCorrelatorGUI.h"


/*=====================
  A class to process and save information about the thermal environment*/
class AnitaNoiseMachine
{
 public:

  /* Constructor */
  AnitaNoiseMachine(const int length = 60);

  const int fifoLength; //one minute of noise averaging

  //do you want to save the interferometric maps?  They are very large.  Also multiple ways to save them
  bool fillMap ;    //save the min bias maps as TH2D (~18kB per event)
  bool fillAvgMap; //save it as the average of fifoLength min bias maps (Won't work if fillMap==true)
  bool fillArray;  //save it as a double array (25% smaller, ~14kB per event)
  
  /* Reset */
  void zeroInternals();

  /* get time avg waveform rms noise value for a channel */
  double getAvgRMSNoise(int phi, int ringi, int poli);

  /* makes an averaged TProfile2D out of the histograms in the mapFifo */
  TProfile2D *getAvgMapNoiseProfile(AnitaPol::AnitaPol_t pol);


  /* updates all the fifos with the current event */
  void updateMachine(UCorrelator::Analyzer *analyzer,FilteredAnitaEvent *filtered);


  /* Moves things into the summary */
  void fillNoiseSummary(AnitaNoiseSummary *noiseSummary); //crabcrabcrab

  /* Fill the mapHistoryVal value in AnitaEventSummary (eventSummary should be mostly filled already) */
  void fillEventSummary(AnitaEventSummary *eventSummary);
  //crab



 private:

  //Makes sure the fifos start at zero
  bool fJustInitialized;

  
  //induvidual update functions which need to all be called at once so the fifos incriment properly
  /* for building up an enormous memory block of histograms from a bunch of events, then making an average 
              also does the double array   */
  void updateAvgMapFifo(UCorrelator::Analyzer *analyzer);
  /* for calculating rms of waveform from a minute average before event capture */
  void updateAvgRMSFifo(FilteredAnitaEvent *filtered);


  //internals for time domain waveform rms fifo
  double *rmsFifo; //where the info is saved
  int rmsFifoPos; //where in the fifo the most recent write was
  bool rmsFifoFillFlag ; //whether you've completely filled the fifo once
  int rmsFifoIndex(int pol,int ringi,int poli,int pos); //so I know where to write in the 1d array

  //internals for interferometric map fifo (probably enormous in memory so maybe make a flag)
  TH2D **mapFifo[NUM_POLS]; //where the info is saved
  int mapFifoPos;  //where in the fifo the most recent write was
  bool mapFifoFillFlag ; //whether you've completely filled the fifo once.
  int mapFifoIndex(); 


  //maybe will be more compressed
  static const int nPhi = 180; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
  static const int nTheta = 100; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
  //where the quickly calculated double array is stored.  should be maxDirections*nPhi*nTheta*NUM_POLS long
  double *rollingMapAvg;
  int rollingMapIndex(int poli,int iPhi,int iTheta); 

  ClassDefNV(AnitaNoiseMachine, 2); 

};
/*------------------*/


#endif
