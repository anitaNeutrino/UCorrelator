#include "Correlator.h" 
#include "AnitaDataset.h" 
#include "FilterStrategy.h" 
#include "BasicFilters.h" 
#include "UsefulAdu5Pat.h" 




void sunCorrelatorExample(int run=352, int event = 60832108,AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal) 
{

  //load the event and run
  AnitaDataset d(run);
  d.getEvent(event); 

  //Filter strategy with just the alfa filter
  FilterStrategy strategy; 
  strategy.addOperation(new ALFAFilter); 

  // make our filtered anita event
  FilteredAnitaEvent  ev(d.useful(), &strategy, d.gps(), d.header()); 

  //since we are not using the non-zoomed histogram, we don't really care about the options here. making it coarser will make that part faster
  UCorrelator::Correlator corr(120, 0,360,40,40,40); 

  //figure out where the sun is 
  double phi_sun,theta_sun; 

  ((UsefulAdu5Pat*) d.gps())->getSunPosition(phi_sun, theta_sun);  //this cast is just because getSunPosition should be const but it's not :( 
  
  //we just want to compute zoomed aorund the sun, but we still have to call compute first otherwise it won't have the event. 
  corr.compute(&ev,pol); 

  //note that the sign convention is reversed between computeZoomed  and the sunPosition (sorry about that!).
  //this is for semi-historical reasons  (the histogram is drawn backwards from the theta sign convention) 
  TH2 * hist = corr.computeZoomed(phi_sun, -theta_sun , 60,0.5,60,0.5); 

  hist->SetTitle("sun centered"); 
  hist->Draw("colz"); 

  //Draw a marker at the sun...
  TMarker * sun_marker = new TMarker(phi_sun,-theta_sun,3); 
  sun_marker->Draw("same"); 



}

