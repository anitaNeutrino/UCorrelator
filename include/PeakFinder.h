#ifndef UCORRELATOR_PEAK_FINDER_H
#define UCORRELATOR_PEAK_FINDER_H

class TH2D; 
#include "AnitaEventSummary.h"

/** \file 
 *  
 *  This file contains various peak finding routines
 *
 */ 
namespace UCorrelator
{

  namespace peakfinder
  {



    /** Stores the location of a rough maximum */
    class RoughMaximum
    {
      public:
        int bin; 
        int val; 
        double x; 
        double y; 
    }; 

    /** Stores the location and uncertainty of a fine maximum */
    class FineMaximum
    {
      public: 
        double x; 
        double y; 
        double sigma_x; 
        double sigma_y;
        double covar; 
        double val; 

        /** copy into pointing hypothesis struct */
        void copyToPointingHypothesis(AnitaEventSummary::PointingHypothesis *p);
    }; 

    /** Finds largest Nmaxima isolated local maxima in a histogram. A maximum is local if it's bigger than any neighboring bins
     * Returns the number found (which might be less than Nmaxima if there are insufficient isolated local maxima) **/ 
    int findIsolatedMaxima(const TH2D* hist, double distance, int Nmaxima, RoughMaximum * maxima, bool use_bin_center = true); 

    /** Legacy interpolation peak finding as Abby did it... basically does interpolation independently in row and column */ 
    void doInterpolationPeakFindingAbby(const TH2D* hist, FineMaximum * peak);

    /** NOT IMPLEMENTED YET. Eventually, this will find the maximum of the bicubic interpolant in the maximum bin */ 
    void doInterpolationPeakFindingBicubic(const TH2D* hist, FineMaximum * peak);

    /** These all do a quadratic fit in the neighbhorhood of the maximum bin. The number denotes the size of the square.
     *  This should be very fast; the matrix decomposition necessary to perform the minimization is statically computed.
     **/ 
    void doPeakFindingQuadratic9(const TH2D* hist, FineMaximum * peak);
    void doPeakFindingQuadratic16(const TH2D* hist, FineMaximum * peak);
    void doPeakFindingQuadratic25(const TH2D* hist, FineMaximum * peak);
    void doPeakFindingQuadratic36(const TH2D* hist, FineMaximum * peak);
    void doPeakFindingQuadratic49(const TH2D* hist, FineMaximum * peak);

    /** 2D gaussian fit to histogram.... this will be pretty slow */ 
    void doPeakFindingGaussian(const TH2D* hist, FineMaximum * peak);

  }

}

#endif
