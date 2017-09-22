#include "SineSubtractCache.h"
#include "FilterStrategy.h"
#include "UCFilters.h"

int main(int argc, char** argv){

  int run = argc >= 2 ? atoi(argv[1]) : 352;
    
  FilterStrategy* strat = new FilterStrategy();
  UCorrelator::fillStrategyWithKey(strat,"sinsub_10_3_ad_2");

  UCorrelator::SineSubtractFilter* ssf = dynamic_cast<UCorrelator::SineSubtractFilter*>(const_cast<FilterOperation*>(strat->getOperation(0)));
  UCorrelator::SineSubtractCache::makeCache(run, ssf);
}

