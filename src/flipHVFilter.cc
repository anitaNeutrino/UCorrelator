#include "FFTtools.h" 
#include "FilterStrategy.h"
#include "flipHVFilter.h" 
#include "BasicFilters.h" 
#include "ResponseManager.h" 
#include "SystemResponse.h" 
#include <map> 
#include "TGraph.h" 
#include "AnalysisWaveform.h" 
#include "RawAnitaHeader.h"
#include "TRandom.h"
#include "TFile.h"
#include "Adu5Pat.h"
#include "DigitalFilter.h"
#include "RawAnitaHeader.h"
#include <math.h>// for isnan
#include "AntennaPositions.h"
#include <ctype.h>
#include <stdio.h>


UCorrelator::flipHVFilter::flipHVFilter() {


}


UCorrelator::flipHVFilter::~flipHVFilter() {


}


void UCorrelator::flipHVFilter::process(FilteredAnitaEvent * ev) {

  processOne(getWf(ev, 44, AnitaPol::kVertical), ev->getHeader(), 44, AnitaPol::kVertical);

}


void UCorrelator::flipHVFilter::processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol) {


}
