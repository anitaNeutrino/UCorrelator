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

  static const int fifoLength = 60; //one minute of noise averaging

  //do you want to save the interferometric maps?  They are very large.  Also multiple ways to save them
  bool fillMap ;    //save the min bias maps as TH2D (~18kB per event)
  bool fillAvgMap; //save it as the average of fifoLength min bias maps (Won't work if fillMap==true)
  bool fillArray ;  //save it as a double array (25% smaller, ~14kB per event)


  /* Constructor */
  AnitaNoiseMachine();
  
  void zeroInternals();

  
  /* for calculating rms of waveform from a minute average before event capture */
  void fillAvgRMSNoise(FilteredAnitaEvent *filtered);
  double getAvgRMSNoise(int phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol);

  /* for building up an enormous memory block of histograms from a bunch of events, then making an average */
  void fillAvgMapNoise(UCorrelator::Analyzer *analyzer);

  /* makes a TProfile2D out of the histograms in the fifo */
  TProfile2D *getAvgMapNoiseProfile(AnitaPol::AnitaPol_t pol);

  /* updates the weird map fifo array that I'm trying to get to work */
  void updateAvgMapNoise();

  // basically just moves things into the summary
  void fillNoiseSummary(AnitaNoiseSummary *noiseSummary); //crabcrabcrab

  //also I added something to eventSummary that I do
  void fillEventSummary(AnitaEventSummary *eventSummary);
  //crab



 private:

  //internals for time domain waveform rms fifo
  double rmsFifo[NUM_PHI][NUM_ANTENNA_RINGS][NUM_POLS][fifoLength]; //where the info is saved
  int rmsFifoPos; //where in the fifo the most recent write was
  bool rmsFifoFillFlag ; //whether you've completely filled the fifo once

  //internals for interferometric map fifo (probably enormous in memory so maybe make a flag)
  TH2D *mapFifo[NUM_POLS][fifoLength]; //where the info is saved
  int mapFifoPos;  //where in the fifo the most recent write was
  bool mapFifoFillFlag ; //whether you've completely filled the fifo once.

  //maybe will be more compressed
  static const int nPhi = 180; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
  static const int nTheta = 100; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
  double avgMaps[NUM_POLS][nPhi][nTheta]; //where the array is stored.  should be nPhi*nTheta*NUM_POLS long

  ClassDefNV(AnitaNoiseMachine, 2); 

};
/*------------------*/


#endif
