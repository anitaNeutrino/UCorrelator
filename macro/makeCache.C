#include "UCFilters.h"
#include "SineSubtractCache.h"

void makeCache(int run = 352){

  FilterStrategy* strat = new FilterStrategy();
  UCorrelator::fillStrategyWithKey(strat,"sinsub_10_3_ad_2");

  auto ssf = dynamic_cast<UCorrelator::SineSubtractFilter*>(const_cast<FilterOperation*>(strat->getOperation(0)));
  UCorrelator::SineSubtractCache::makeCache(run, ssf);
}
