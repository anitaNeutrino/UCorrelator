#include "Baseline.h" 
#include "TGraph.h" 
#include <cstdio>
#include "TSystem.h" 
#include "AnitaDataset.h"
#include "UsefulAnitaEvent.h"
#include "FFTtools.h"
#include "RawAnitaHeader.h"
#include "Flags.h"
#include "TFile.h" 
#include "FFTWComplex.h" 

#define NSAMPLES 256  
#define NOMINAL_DT 1./2.6 


static int makeBaselines(int run, TGraph ** hpol, TGraph ** vpol, int N = 5000) 
{
  char *  datadir = getenv("ANITA_ROOT_DATA"); 
  if (!datadir) 
  {
    fprintf(stderr,"The environmental variable ANITA_ROOT_DATA must be defined! Aborting."); 
    exit(1); 
  }


  AnitaDataset d(run); 

  int nevents = 0; 
  int i = 10; 


  while (nevents < N && i < d.N()) 
  {
    d.getEntry(i); 

    // skip non-RF triggers and a small fraction of times? is this really right? This is what Abby has though
    if (d.header()-> trigType != 1 && d.header()->triggerTimeNs <= 1e6)  
       continue; 

    //skip saturated events 
    if (UCorrelator::flags::checkSaturation(d.useful())) 
      continue; 

    for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
    {
       if (nevents == 0)
       {
           hpol[ant] = new TGraph(NSAMPLES/2+1); 
           vpol[ant] = new TGraph(NSAMPLES/2+1); 
       }
       
       TGraph * gh = d.useful()->getGraph(ant, AnitaPol::kHorizontal); 
       TGraph * gv = d.useful()->getGraph(ant, AnitaPol::kVertical); 

       TGraph * igh = FFTtools::getInterpolatedGraph(gh, NOMINAL_DT); 
       TGraph * igv = FFTtools::getInterpolatedGraph(gv, NOMINAL_DT); 

       igh->Set(NSAMPLES); 
       igv->Set(NSAMPLES); 

       FFTWComplex * fft_h = FFTtools::doFFT(igh->GetN(), igh->GetY()); 
       FFTWComplex * fft_v = FFTtools::doFFT(igv->GetN(), igv->GetY()); 

       for (int j = 0; j < NSAMPLES/2+1; j++) 
       {
         if (nevents == 0) 
         {
           hpol[ant]->GetX()[j] = j  / (NOMINAL_DT * NSAMPLES);  //GHz
           vpol[ant]->GetX()[j] = j  / (NOMINAL_DT * NSAMPLES);  //GHz
         }

         hpol[ant]->GetY()[j] += fft_h[i].getAbs() / (NSAMPLES/2.+1)  * (j == 0 || j == NSAMPLES/2 ? 1 : 2); 
         vpol[ant]->GetY()[j] += fft_v[i].getAbs() / (NSAMPLES/2.+1)  * (j == 0 || j == NSAMPLES/2 ? 1 : 2); 

       }

       delete [] fft_h; 
       delete [] fft_v; 
       delete gh; 
       delete gv; 
       delete igh; 
       delete igv; 

    }

    nevents++; 
  }

  for (int ant = 0; ant < NUM_SEAVEYS; ant++)
  {
    for (int i = 0; i < NSAMPLES / 2 + 1; i++) 
    {
      hpol[ant]->GetY()[i] /= nevents; 
      vpol[ant]->GetY()[i] /= nevents; 
    }
  }


  return nevents; 
}


void UCorrelator::Baseline::saveToDir(const char * dir) 
{

  gSystem->mkdir(dir,true); 
  TFile f(TString::Format("%s/%d_%d.root", dir, run,navg),"RECREATE"); 
  for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
  {
    hpol[ant]->Write(TString::Format("h%d",ant)); 
    vpol[ant]->Write(TString::Format("v%d",ant)); 
  }
}

UCorrelator::Baseline::Baseline(int run, int navg, const char * persistdir) 
  : navg(navg) , run(run)
{
  if (persistdir)  
  {
    TFile f(TString::Format("%s/%d_%d.root", persistdir, run,navg)); 
    if (f.IsOpen())
    {
      for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
      {
        TGraph * found_hpol = (TGraph*) f.Get(TString::Format("h%d",ant)); 
        hpol[ant] = new TGraph(*found_hpol);  

        TGraph * found_vpol = (TGraph*) f.Get(TString::Format("v%d",ant)); 
        vpol[ant] = new TGraph(*found_vpol);  
      }
      return; 
    }
  }

  makeBaselines(run, hpol, vpol, navg); 


  hpol_avg = new TGraph(NSAMPLES/2+1); 
  vpol_avg = new TGraph(NSAMPLES/2+1); 

  memcpy(hpol_avg->GetX(), hpol[0]->GetX(), sizeof(double) * hpol_avg->GetN());
  memcpy(vpol_avg->GetX(), vpol[0]->GetX(), sizeof(double) * vpol_avg->GetN());

  for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
  {
    for (int i =0; i < NSAMPLES/2 +1; i++) 
    {
      hpol_avg->GetY()[i] += hpol[ant]->GetY()[i] / NUM_SEAVEYS; 
      vpol_avg->GetY()[i] += vpol[ant]->GetY()[i] / NUM_SEAVEYS; 
    }
  }

  if (persistdir) saveToDir(persistdir); 
}
            
