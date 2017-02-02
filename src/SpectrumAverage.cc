

#include "SpectrumAverage.h" 
#include "TFile.h" 
#include "TSystem.h" 
#include "AnitaVersion.h" 
#include "AnalysisWaveform.h" 
#include "FilterStrategy.h"
#include "TH2.h" 
#include "RawAnitaHeader.h" 
#include "AnitaDataset.h" 
#include "TCut.h" 
#include "FilteredAnitaEvent.h" 
#include "BasicFilters.h"



static ALFAFilter alfa; 




static int computeAverage(int run, int nsecs, const char * selection, double max_r, TH2 * avgs[NUM_SEAVEYS][2]) 
{

  AnitaDataset d(run); 
  if (selection) d.setCut(selection); 


  //figure out how long it is. 

  selection ? d.firstInCut() : d.first(); 
  double startTime = d.header()->triggerTime; 
  selection ? d.lastInCut() : d.last(); 
  double endTime = d.header()->triggerTime+1; 

  int nbins = (startTime - endTime) / nsecs; 


  TH1I norm("specnorm","Normalization",nbins,startTime,endTime); 

  for (int ant = 0; ant < NUM_SEAVEYS; ant++)
  {
    for (int pol = 0; pol < 2; pol++) 
    {

      avgs[ant][pol] = new TH2D(TString::Format("specavg_%d_%d_%d_%d",run,nsecs,ant,pol), 
                         TString::Format("Spectrum Average ant=%d, pol=%d, run %d, nsec=%d, sel=%s", ant,pol,run,nsecs,selection),
                         nbins,startTime,endTime, 
                         131,0,1.31); //TODO don't hardcode here
    }
  }



                       
  FilterStrategy str; 
  if (AnitaVersion::get() == 3) str.addOperation(&alfa); 


  int N = selection ? d.NInCut():  d.N(); 

  for (int i = 0; i < N; i++)
  {
    selection ? d.nthInCut(i) : d.getEntry(i);
    FilteredAnitaEvent ev(d.useful(), &str, d.gps(), d.header()); 
    
    //cut out possible blasts 
    double max_ratio_hpol, max_ratio_vpol; 
    ev.getMinMaxRatio(AnitaPol::kHorizontal, &max_ratio_hpol,0,0,0);
    if (max_ratio_hpol > max_r) continue; 
    ev.getMinMaxRatio(AnitaPol::kVertical, &max_ratio_vpol,0,0,0); 
    if (max_ratio_vpol > max_r) continue; 


    double t= d.header()->triggerTime; 
    norm.Fill(t); 
    int xbin = norm.FindFixBin(t); 

#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for 
#endif
    for (int j = 0; j < NUM_SEAVEYS * 2; j++) 
    {
      int ant = j /2; 
      int pol = j %2; 
      const AnalysisWaveform * wf = ev.getFilteredGraph(ant,AnitaPol::AnitaPol_t(pol)); 

      const TGraphAligned * g = wf->power(); 
      for (int k = 0; k<avgs[ant][pol]->GetNbinsY() ; k++)
      {
        avgs[ant][pol]->SetBinContent(xbin,k, avgs[ant][pol]->GetBinContent(xbin,k) + g->GetY()[k]); 
      }
    }
  }

  return 0; 
}



void UCorrelator::SpectrumAverage::saveToDir(const char * dir) 
{

  gSystem->mkdir(dir,true); 
  TFile f(TString::Format("%s/%d_%d.root", dir, run,nsecs),"RECREATE"); 

  for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
  {
   avgs[ant][0]->Write(TString::Format("h%d",ant)); 
   avgs[ant][1]->Write(TString::Format("v%d",ant)); 
  }
}

UCorrelator::SpectrumAverage::SpectrumAverage(int run, int nsecs, const char * persistdir,
     const char * selection, double max_bottom_top_ratio) 
  : nsecs(nsecs) , run(run)
{
  bool foundit = false;

  const char * check_dir = persistdir ? persistdir 
                                      :  getenv("UCORRELATOR_SPECAVG_DIR") ?: ""; 


  if (check_dir)
  {

    TFile f(TString::Format("%s/%d_%d.root", check_dir, run,nsecs)); 
    if (f.IsOpen())
    {
      for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
      {
        TH2D * found_hpol = (TH2D*) f.Get(TString::Format("h%d",ant)); 
        avgs[ant][0] = new TH2D(*found_hpol);  

        TH2D * found_vpol = (TH2D*) f.Get(TString::Format("v%d",ant)); 
        avgs[ant][1] = new TH2D(*found_vpol);  
      }
      foundit = true;
    }
  }

  if (!foundit) 
  {
    printf("Didn't find spectrum_averages... creating!\n"); 
    if (!persistdir) 
    {
      printf("Define UCORRELATOR_SPECAVG_DIR to persist somewhere.\n"); 
    }
    computeAverage(run, nsecs, selection, max_bottom_top_ratio, avgs); 
    if (persistdir) saveToDir(persistdir); 
  }

}
            


