#ifndef UCORRELATOR_IMAGE_TOOLS_H
#define UCORRELATOR_IMAGE_TOOLS_H

/** Image tools for working with TH2*'s
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 **/


class TH1; 
class TH2; 


namespace UCorrelator
{

  /*=====
    TH1::GetRMS(3) stopped working so here is a function to just calculate it by hand
   */
  double getZRMS(const TH2*);



  namespace image
  {

    TH1 * getPctileProjection(const TH2 * h, int axis = 1, double pct = 0.5, bool ignoreEmpty = true); 


    enum InterpolationType
    {
      NEAREST, 
      BILINEAR, 
      BICUBIC
    };


    enum  InterpolationEdgeBehavior 
    {
      OVERFLOW_BIN,
      EXTEND,
      PERIODIC
    }; 

 
    /**
     * Interpolate histogram h at values x and y 
     * 
     * Currently, we support nearest neighbor, bilinear, and bicubic interpolation.
     * Note however that bicubic only works for histograms with equal bin sizes (if the bin sizes are not equal, we delegate to bilinear).
     *
     * The edge behavior can be set to either overflow bin (likely 0), extend (to extend the non-overflow bin value infinitely) or periodic for each axis. 
     * 
     * You can also toggle whether or not the histogram value represents the center (default) or the lower edge. 
     *
     * @param h The histogram to interpolate
     * @param x The x value at which to interpolate
     * @param y The y value at which to interpolate
     * @param type The interpolation type 
     * @param edge_x x-coord boundary behavior
     * @param edge_y y-coord boundary behavior
     * @param use_bin_centers true to use bin centers
     *
     *
     * */ 
    double interpolate(const TH2 * h, double x, double y, InterpolationType type = BILINEAR, InterpolationEdgeBehavior edge_x= PERIODIC, InterpolationEdgeBehavior edge_y = OVERFLOW_BIN, bool use_bin_centers = true); 

  }; 
}; 

#endif
