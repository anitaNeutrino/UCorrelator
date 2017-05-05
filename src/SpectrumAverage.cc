#include "SpectrumAverage.h" 
#include "AnitaVersion.h" 
#include "TFile.h" 
#include "UCImageTools.h" 
#include "TSystem.h" 
#include "AnitaVersion.h" 
#include "AnalysisWaveform.h" 
#include "FilterStrategy.h"
#include "TH2.h" 
#include "RawAnitaHeader.h" 
#include "AnitaDataset.h" 
#include "TCut.h" 
#include "FilteredAnitaEvent.h" 
#include "TMutex.h" 
#include "BasicFilters.h"
#include <math.h>// for isnan

static ALFAFilter alfa; 


//#define DEBUG_SPEC_AVG 


UCorrelator::SpectrumAverage::~SpectrumAverage()
{

  for (int i = 0; i < NUM_SEAVEYS; i++) 
  {
    for (int j = 0; j < 2; j++) 
    {
      delete avgs[i][j]; 
      if (peakiness[i][j]) delete peakiness[i][j]; 

    }
  }

}

int UCorrelator::SpectrumAverage::computeAverage(const char * selection, double max_r) 
{

  AnitaDataset d(run); 
#ifdef MULTIVERSION_ANITA_ENABLED
  if (selection) d.setCut(selection); 
#else
  if (selection) fprintf(stderr,"selection not supported with this version of AnitaDataset\n"); 
#endif



  //figure out how long it is. 

#ifdef MULTIVERSION_ANITA_ENABLED
  selection ? d.firstInCut() : d.first(); 
#else
  d.getEntry(0); 
#endif
  double startTime = d.header()->triggerTime; 
#ifdef MULTIVERSION_ANITA_ENABLED
  selection ? d.lastInCut() : d.last(); 
#else
  d.getEntry(d.N()-1);
#endif
  double endTime = d.header()->triggerTime+1; 

  int nbins = (endTime-startTime) / nsecs; 

//  printf("%d\n", nbins); 

  TH1I norm("specnorm","Normalization",nbins,startTime,endTime); 

  for (int ant = 0; ant < NUM_SEAVEYS; ant++)
  {
    for (int pol = 0; pol < 2; pol++) 
    {
      TString name; 
      name.Form("specavg_%d_%d_%d_%d",run,nsecs,ant,pol); 
      TString title; 
      title.Form("Spectrum Average ant=%d, pol=%d, run %d, nsec=%d, sel=%s", ant,pol,run,nsecs,selection); 

        
      avgs[ant][pol] = new TH2D(name, title, 
                        131,0,1.31,//TODO don't hardcode here
                        nbins, startTime,endTime) ; 
                         

      avgs[ant][pol]->SetDirectory(0); 
      avgs[ant][pol]->GetYaxis()->SetTitle("Time"); 
      avgs[ant][pol]->GetYaxis()->SetTimeDisplay(true); 
      avgs[ant][pol]->GetXaxis()->SetTitle("Frequency"); 
    }
  }



                       
  FilterStrategy str; 
  if (AnitaVersion::get() == 3) str.addOperation(&alfa); 


#ifdef MULTIVERSION_ANITA_ENABLED
  int N = selection ? d.NInCut():  d.N(); 
#else
  int N = d.N(); 
#endif

#ifdef DEBUG_SPEC_AVG
  N = d.N()/100; 
#endif
  for (int i = 0; i < N; i++)
  {
#ifdef MULTIVERSION_ANITA_ENABLED
    selection ? d.nthInCut(i) : d.getEntry(i);
#else
    d.getEntry(i); 
#endif 
    FilteredAnitaEvent ev(d.useful(), &str, d.gps(), d.header()); 
    
    //cut out possible blasts 
    double max_ratio_hpol, max_ratio_vpol; 
    ev.getMinMaxRatio(AnitaPol::kHorizontal, &max_ratio_hpol,0,0,0);
    if (max_ratio_hpol > max_r) continue; 
    ev.getMinMaxRatio(AnitaPol::kVertical, &max_ratio_vpol,0,0,0); 
    if (max_ratio_vpol > max_r) continue; 


    double t= d.header()->triggerTime; 
    norm.Fill(t); 
    int ybin = norm.FindFixBin(t); 

//#ifdef UCORRELATOR_OPENMP
//#pragma omp parallel for 
//#endif
    for (int j = 0; j < NUM_SEAVEYS * 2; j++) 
    {
      int ant = j /2; 
      int pol = j %2; 
      const AnalysisWaveform * wf = ev.getFilteredGraph(ant,AnitaPol::AnitaPol_t(pol)); 

      const TGraphAligned * g = wf->power(); 
      for (int k = 0; k<avgs[ant][pol]->GetNbinsX() ; k++)
      {
        avgs[ant][pol]->SetBinContent(k,ybin, avgs[ant][pol]->GetBinContent(k,ybin) + g->GetY()[k]); 
      }
    }

#ifdef DEBUG_SPEC_AVG
    if (i % 50 == 0) 
      printf("\n%d",i); 
    else
      printf("."); 
#endif
  }

  for (int j = 0; j < NUM_SEAVEYS * 2; j++) 
  {
    int ant = j /2; 
    int pol = j %2; 
    for (int ybin = 1; ybin <= avgs[ant][pol]->GetNbinsY(); ybin++)
    {
      for (int k = 0; k<avgs[ant][pol]->GetNbinsX() ; k++)
      {
        if (norm.GetBinContent(ybin))
          avgs[ant][pol]->SetBinContent(k,ybin, avgs[ant][pol]->GetBinContent(k,ybin) / norm.GetBinContent(ybin)); 
      }
    }
  }

  return 0; 
}



void UCorrelator::SpectrumAverage::saveToDir(const char * dir) 
{

  gSystem->mkdir(dir,true); 

  TString fstr; 
  fstr.Form("%s/%d_%d.root", dir, run,nsecs); 
  TFile f(fstr.Data(),"RECREATE"); 
  f.cd(); 

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
                                      :  getenv("UCORRELATOR_SPECAVG_DIR") ? : 0; 


  if (check_dir)
  {

    TFile f(TString::Format("%s/%d_%d.root", check_dir, run,nsecs)); 
    if (f.IsOpen())
    {
      for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
      {
        TH2D * found_hpol = (TH2D*) f.Get(TString::Format("h%d",ant)); 
        avgs[ant][0] = new TH2D(*found_hpol);  
        avgs[ant][0]->SetDirectory(0); 

        TH2D * found_vpol = (TH2D*) f.Get(TString::Format("v%d",ant)); 
        avgs[ant][1] = new TH2D(*found_vpol);  
        avgs[ant][1]->SetDirectory(0); 
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
    computeAverage( selection, max_bottom_top_ratio); 
    if (persistdir) saveToDir(persistdir); 
  }

  memset(peakiness,0,sizeof(peakiness)); 
}
            

TH1* UCorrelator::SpectrumAverage::getSpectrumAverage(AnitaPol::AnitaPol_t pol, int ant, double t, bool db)  const
{

  //figure out which bin we are in 
   int bin =  avgs[ant][pol]->GetYaxis()->FindFixBin(t); 

   if (bin == 0 || bin == avgs[ant][pol]->GetNbinsY()+1) return 0 ; 

   TString title; 
   title.Form("Ant%d %c-Pol Average at t=%g", ant, pol == AnitaPol::kHorizontal ? 'H' : 'V', t); 
   TString name;
   name.Form("specavg_%d_%d_%g",ant,pol,t); 
   TH1D * answer =  new TH1D(name,title, 131,0,1.31); //todo don't hardcode 

   for (int i = 1; i <= answer->GetNbinsX(); i++)
   {
     answer->SetBinContent(i, avgs[ant][pol]->GetBinContent(i,bin)); 
   }

   answer->GetXaxis()->SetTitle("Frequency"); 
   answer->GetYaxis()->SetTitle(db ? "Power (dBish)" :"Power (linear)"); 


   if (db) 
   {
     for (int i = 1; i <= answer->GetNbinsX(); i++) 
     {
       answer->SetBinContent(i, 10 * TMath::Log10(answer->GetBinContent(i))); 
     }
   }


   return answer;
  

}

TH1 *UCorrelator::SpectrumAverage::getSpectrumPercentile(AnitaPol::AnitaPol_t pol, int ant, double pct , bool db ) const
{

  TH1 * answer = UCorrelator::image::getPctileProjection( avgs[ant][pol], 1, pct); 

  answer->GetXaxis()->SetTitle("Frequency"); 
  answer->GetYaxis()->SetTitle(db ? "Pctile Power (dBish)" :"Pctile Power (linear)"); 

   answer->SetTitle(TString::Format("Ant%d %c-Pol %gth Pctile", ant, pol == AnitaPol::kHorizontal ? 'H' : 'V', pct*100)); 

   if (db) 
   {
     for (int i = 1; i <= answer->GetNbinsX(); i++) 
     {
       answer->SetBinContent(i, 10 * TMath::Log10(answer->GetBinContent(i))); 
     }
   }


   return answer;
 
}


double UCorrelator::SpectrumAverage::getStartTime() const 
{

  return avgs[0][0] ? avgs[0][0]->GetYaxis()->GetXmin() : -1; 
}

double UCorrelator::SpectrumAverage::getEndTime() const 
{
  return avgs[0][0] ? avgs[0][0]->GetYaxis()->GetXmax() : -1;   
}


static const UCorrelator::SpectrumAverage* defaultThermalAvg = 0; 
const UCorrelator::SpectrumAverage* UCorrelator::SpectrumAverage::defaultThermal() 
{
  if (defaultThermalAvg) return defaultThermalAvg; 

  //TODO make this a4 compatible as well 
   if (AnitaVersion::get() == 4) 
   {
     fprintf(stderr,"warning: using default terminated thermal spectrum for A3 for peakiness\n"); 
   }

   TString dir; 
   dir.Form("%s/share/UCorrelator/terminated_noise/", getenv("ANITA_UTIL_INSTALL_DIR")); 

   defaultThermalAvg = new SpectrumAverage(11382,60, dir.Data()); 
   return defaultThermalAvg; 
}



void UCorrelator::SpectrumAverage::computePeakiness(const SpectrumAverage * thermalSpec, double fractionForNormalization) 
{

  if (!thermalSpec) thermalSpec = defaultThermal(); 

  for (int ant  = 0; ant < NUM_SEAVEYS; ant++)
  {
//    printf("%d\n",ant);
    for (int ipol = 0; ipol < 2; ipol++) 
    {

//      printf("%d %d\n",ant,ipol);

      //get median spectrum 
      TH1 * median = getSpectrumPercentile(AnitaPol::AnitaPol_t(ipol), ant, 0.5); 
      TH1 * thermal = thermalSpec->getSpectrumPercentile(AnitaPol::AnitaPol_t(ipol),ant,0.5);

      // now we have to figure out scale. To do this, we compute the mean of the middle frac of points. 
      int index_spec[median->GetNbinsX()]; 
      int index_therm[thermal->GetNbinsX()]; 

      TMath::Sort(thermal->GetNbinsX(),((TH1D*) thermal)->GetArray()+1, index_therm); 
      TMath::Sort(median->GetNbinsX(), ((TH1D*)median)->GetArray()+1, index_spec); 

      double sum_spec = 0; 
      double sum_therm =0;  

      for (int i = 0; i <= thermal->GetNbinsX()*fractionForNormalization; i++) 
      {
        sum_therm+=thermal->GetBinContent(index_therm[ (int) (i + (thermal->GetNbinsX() * (0.5-fractionForNormalization)))  ]); 
      }

      for (int i = 0; i <= median->GetNbinsX()*fractionForNormalization; i++) 
      {
        sum_spec+=median->GetBinContent(index_spec[ (int) (i + median->GetNbinsX() * (0.5-fractionForNormalization))]); 
      }

      sum_therm /= thermal->GetNbinsX(); 
      sum_spec /= median->GetNbinsX(); 
      double ratio = sum_spec / sum_therm; 

      TString name; 
      name.Form("peakiness_%d_%d\n", ant,ipol); 
      TString title; 
      title.Form("peakiness ant=%d pol=%d\n", ant,ipol); 

      peakiness[ant][ipol] = new TH2D(name,title,
                                       avgs[ant][ipol]->GetNbinsX(), avgs[ant][ipol]->GetXaxis()->GetXmin(), avgs[ant][ipol]->GetXaxis()->GetXmax(), 
                                       avgs[ant][ipol]->GetNbinsY(), avgs[ant][ipol]->GetYaxis()->GetXmin(), avgs[ant][ipol]->GetYaxis()->GetXmax());  


      peakiness[ant][ipol]->SetDirectory(0); 

      for (int ii = 1; ii < peakiness[ant][ipol]->GetNbinsX(); ii++)
      {
        for (int jj = 1; jj < peakiness[ant][ipol]->GetNbinsY(); jj++)
        {
          peakiness[ant][ipol]->SetBinContent(ii,jj, avgs[ant][ipol]->GetBinContent(ii,jj)  / ( ratio * thermal->GetBinContent(ii))); 
          if (std::isnan(peakiness[ant][ipol]->GetBinContent(ii,jj)))
              peakiness[ant][ipol]->SetBinContent(ii,jj); 
        }
      }

      delete median; 
      delete thermal; 

    }
  }

}

UCorrelator:: SpectrumAverageLoader::SpectrumAverageLoader(const char * the_dir, int secs) 
  : dir(the_dir), nsecs(secs)
{
  spec = 0; 
}

static TMutex m; 

const UCorrelator::SpectrumAverage* UCorrelator::SpectrumAverageLoader::avg(double t) const
{

  if (spec && t >= spec->getStartTime() && t <= spec->getEndTime() ) return spec; 

  m.Lock(); 

  //double check 
  if (spec && t >= spec->getStartTime() && t <= spec->getEndTime() ) 
  {
    m.UnLock(); 
    return spec; 
  }

  int run = AnitaDataset::getRunAtTime(t); 
  if (spec) delete spec; 
  spec = new SpectrumAverage(run,nsecs, dir); 
  spec->computePeakiness(); 
  m.UnLock(); 
  
  return spec; 

}
