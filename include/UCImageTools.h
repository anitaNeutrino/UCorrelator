#ifndef UCORRELATOR_IMAGE_TOOLS_H
#define UCORRELATOR_IMAGE_TOOLS_H

/** Image tools for working with TH2*'s
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 *  Also other things!  Ben Rotter <BenJRotter@gmail.com>
 **/


class TH1; 
class TH2; 
class TProfile2D;

namespace UCorrelator
{

  /*=====
    TH1::GetRMS(3) stopped working so here is a function to just calculate it by hand
    If lims is passed, the rms will only be computed within bins  [ x0, x1, y0, y1] 
   */
  double getZRMS(const TH2*, const int * lims = 0);

  /*====
    Need a way to easily rotate the maps to account for heading.  Makes a new map since the input is probably const
    Always returns something that goes from 0->360 with 0/360 being north
  */
  TH2* rotateHistogram(const TH2* inHist,double rotate);

  /*===
    Tools for averaging together interferometric maps, or any histograms really:
    TH2toTProfile2D - "upgrades" a TH2 to a TProfile2D, which can do more stuff. Doesn't copy any of the data, just the frame.
   */
  TProfile2D *TH2toTProfile2D(TH2* inTH2);

  /*===
    Tools for averaging together interferometric maps, or any histograms really:
    fillTProfile2DWithTH2 - "overlays" a TH2 onto a TProfile2D, which will average them together
   */
  void fillTProfile2DWithTH2(TProfile2D *prof, TH2* hist);


  namespace image
  {

    TH1 * getPctileProjection(const TH2 * h, int axis = 1, double pct = 0.5, bool ignoreEmpty = true, const TH1 * norm = 0); 

    /** Converts a histogram to percent below level */ 
    TH1 * makePctileHist(const TH1 * h, const char * name = "_pct", bool invert = false,  int npctilebins = 1000, bool integral=true); 
    
    enum AreaMode 
    {
      AreaMode_Contour, 
      AreaMode_BinCartesian,
      AreaMode_BinSpherical
    }; 
//    double getAreaAbove(const TH2 *h, double val, int mode) 

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
