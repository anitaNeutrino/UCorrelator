#include "ShapeParameters.h" 
#include "TGraph.h" 


double UCorrelator::shape::getRiseTime(const TGraph * g, double min, double max) 
{
  int ifirst = -1; 
  int ilast = -1; 


  for (int i = 0; i < g->GetN(); i++) 
  {

    if (g->GetY()[i] >= min && ifirst < 0)
    {
      ifirst = i; 
    }
    if (g->GetY()[i] >= max && ilast < 0)
    {
      ilast = i; 
      break; 
    }
  }

  if (ifirst < 0 || ilast < 0) return -1; 

  return g->GetX()[ilast] - g->GetX()[ifirst]; 
}

double UCorrelator::shape::getFallTime(const TGraph * g, double min, double max) 
{
  int ifirst = -1; 
  int ilast = -1; 


  for (int i = g->GetN()-1; i >= 0; i--)
  {

    if (g->GetY()[i] >= min && ilast < 0)
    {
      ilast = i; 
    }
    if (g->GetY()[i] >= max && ifirst < 0)
    {
      ifirst = i; 
      break; 
    }
  }

  if (ifirst < 0 || ilast < 0) return -1; 

  return g->GetX()[ilast] - g->GetX()[ifirst]; 
}


double UCorrelator::shape::getWidth(const TGraph * g, double val, int * start, int * end)
{

  int ifirst = -1; 
  int ilast = -1; 

  for (int i = 0; i < g->GetN(); i++) 
  {
    if (g->GetY()[i] >=val) 
    {
      ifirst = i; 
      break; 
    }
  }
  for (int i = g->GetN()-1; i >= 0; i--)
  {
    if (g->GetY()[i] >=val) 
    {
      ilast = i; 
      break; 
    }
  }

  if (ifirst < 0 || ilast< 0) return -1; 

  double t0 = g->GetX()[ifirst]; 
  double t1 = g->GetX()[ilast]; 

  if (start) *start = ifirst; 
  if (end) *end = ilast; 

  return t1-t0; 

}
