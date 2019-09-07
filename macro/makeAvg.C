#include "AnitaVersion.h" 
void makeAvg(int run, int nsecs=10, int anita = 3)
{
  UCorrelator::TimeDependentAverage avg(run, nsecs, "timeavg"); 
}
