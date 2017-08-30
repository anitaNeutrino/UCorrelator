#include "TimeDependentAverage.h" 
#include "AnitaVersion.h" 
#include "TFile.h" 
#include "UCImageTools.h" 
#include "TSystem.h"
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


#define DEBUG_SPEC_AVG 


UCorrelator::TimeDependentAverage::~TimeDependentAverage()
{

  for (int i = 0; i < NUM_SEAVEYS; i++) 
  {
    for (int j = 0; j < 2; j++) 
    {
      delete avgs[i][j]; 
      delete avgs_minbias[i][j]; 
      delete rms[i][j]; 
      if (peakiness[i][j]) delete peakiness[i][j]; 
      if (peakiness_minbias[i][j]) delete peakiness_minbias[i][j]; 

    }
  }

  delete nblasts; 
  delete norms; 
  delete norms_minbias; 
}

int UCorrelator::TimeDependentAverage::computeAverage(double max_r, int min_norm, double max_power) 
{
  AnitaDataset d(run); 




  //figure out how long it is. 

  d.getEntry(0); 
  double startTime = d.header()->triggerTime; 
  d.getEntry(d.N()-1);
  double endTime = d.header()->triggerTime+1; 

  int nbins = (endTime-startTime) / nsecs; 

  TString name; 
  TString title; 
  name.Form("nblasts%d_%d",run,nsecs); 
  title.Form("N Blast Candidates run %d, nsec=%d", run,nsecs); 
  nblasts = new TH1I(name, title, nbins, startTime,endTime) ; 
  nblasts->SetDirectory(0); 
  nblasts->GetYaxis()->SetTitle("NBlasts"); 
  nblasts->GetXaxis()->SetTitle("Time"); 
  nblasts->GetXaxis()->SetTimeDisplay(true); 
  nblasts->GetXaxis()->SetTimeOffset(0); 

  name.Form("norm_%d_%d",run,nsecs); 
  title.Form("RF Normalization run %d, nsec=%d", run,nsecs); 
  
  norms = new TH1I(name, title,nbins,startTime,endTime); 
  norms->SetDirectory(0); 
  norms->GetYaxis()->SetTitle("NEvents"); 
  norms->GetXaxis()->SetTitle("Time"); 
  norms->GetXaxis()->SetTimeDisplay(true); 
  norms->GetXaxis()->SetTimeOffset(0); 

  name.Form("norm_minbias_%d_%d",run,nsecs); 
  title.Form("MinBias Normalization run %d, nsec=%d", run,nsecs); 
  
  norms_minbias = new TH1I(name, title,nbins,startTime,endTime); 
  norms_minbias->SetDirectory(0); 
  norms_minbias->GetYaxis()->SetTitle("NEvents"); 
  norms_minbias->GetXaxis()->SetTitle("Time"); 
  norms_minbias->GetXaxis()->SetTimeDisplay(true); 
  norms_minbias->GetXaxis()->SetTimeOffset(0); 


//  printf("%d\n", nbins); 

  for (int ant = 0; ant < NUM_SEAVEYS; ant++)
  {
    for (int pol = 0; pol < 2; pol++) 
    {

      name.Form("specavg_%d_%d_%d_%d",run,nsecs,ant,pol); 
      title.Form("RF Spectrum Average ant=%d, pol=%d, run %d, nsec=%d", ant,pol,run,nsecs); 

        
      avgs[ant][pol] = new TH2F(name, title, 
                        131,0,1.31,//TODO don't hardcode here
                        nbins, startTime,endTime) ; 
      avgs[ant][pol]->SetDirectory(0); 
      avgs[ant][pol]->GetYaxis()->SetTitle("Time"); 
      avgs[ant][pol]->GetYaxis()->SetTimeDisplay(true); 
      avgs[ant][pol]->GetYaxis()->SetTimeOffset(0); 
      avgs[ant][pol]->GetXaxis()->SetTitle("Frequency"); 


      name.Form("specavg_minbias_%d_%d_%d_%d",run,nsecs,ant,pol); 
      title.Form("MinBias Spectrum Average ant=%d, pol=%d, run %d, nsec=%d", ant,pol,run,nsecs); 

      avgs_minbias[ant][pol] = new TH2F(name, title, 
                        131,0,1.31,//TODO don't hardcode here
                        nbins, startTime,endTime) ; 
      avgs_minbias[ant][pol]->SetDirectory(0); 
      avgs_minbias[ant][pol]->GetYaxis()->SetTitle("Time"); 
      avgs_minbias[ant][pol]->GetYaxis()->SetTimeDisplay(true); 
      avgs_minbias[ant][pol]->GetYaxis()->SetTimeOffset(0); 
      avgs_minbias[ant][pol]->GetXaxis()->SetTitle("Frequency"); 

      name.Form("rms_%d_%d_%d_%d",run,nsecs,ant,pol); 
      title.Form("MinBias RMS ant=%d, pol=%d, run %d, nsec=%d", ant,pol,run,nsecs); 
      rms[ant][pol] = new TH1D(name, title, nbins, startTime,endTime) ; 
      rms[ant][pol]->SetDirectory(0); 
      rms[ant][pol]->GetYaxis()->SetTitle("RMS (mVish)"); 
      rms[ant][pol]->GetXaxis()->SetTitle("Time"); 
      rms[ant][pol]->GetXaxis()->SetTimeDisplay(true); 
      rms[ant][pol]->GetXaxis()->SetTimeOffset(0); 

 


    }
  }



                       
  FilterStrategy str; 
  if (AnitaVersion::get() == 3) str.addOperation(&alfa); 


  int N = d.N(); 

   N = d.N(); //For all the events.
  for (int i = 0; i < N; i++)
  {
    d.getEntry(i); 
    FilteredAnitaEvent ev(d.useful(), &str, d.gps(), d.header()); 
    
    double t= d.header()->triggerTime; 

    //cut out possible blasts. This is perhaps a bit aggressive, but it's worth it 
    double max_ratio_hpol, max_ratio_vpol; 
    ev.getMinMaxRatio(AnitaPol::kHorizontal, &max_ratio_hpol,0,0,0);
    if (max_ratio_hpol > max_r || ev.getAveragePower(AnitaPol::kHorizontal)  > max_power) 
    {
      nblasts->Fill(t); 
      continue; 
    }
    ev.getMinMaxRatio(AnitaPol::kVertical, &max_ratio_vpol,0,0,0); 
    if (max_ratio_vpol > max_r || ev.getAveragePower(AnitaPol::kVertical) > max_power)
    {
      nblasts->Fill(t); 
      continue; 
    }

    bool isRF = d.header()->trigType & 1; 

    (isRF ? norms : norms_minbias)->Fill(t); 
    int tbin = norms->FindFixBin(t); 

#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for 
#endif
    for (int j = 0; j < NUM_SEAVEYS * 2; j++) 
    {
      int ant = j /2; 
      int pol = j %2; 
      const AnalysisWaveform * wf = ev.getFilteredGraph(ant,AnitaPol::AnitaPol_t(pol)); 


      if (!isRF) 
      {
        rms[ant][pol]->SetBinContent(tbin, rms[ant][pol]->GetBinContent(tbin) + wf->even()->GetRMS(2)); 
      }

      const TGraphAligned * g = wf->power(); 
      for (int k = 0; k<avgs[ant][pol]->GetNbinsX() ; k++)
      {
        TH2F * h = isRF ? avgs[ant][pol] : avgs_minbias[ant][pol]; 
        h->SetBinContent(k,tbin,h->GetBinContent(k,tbin) + g->GetY()[k]); 
      }
    }

    if (i % 50 == 0)
    {
      printf("\r%d/%d",i,N);
    }
    else {
      printf(".");
      fflush(stdout);
    }
    
  }
  printf("..Done\n"); 


  for (int isRF = 0; isRF <=1; isRF++)
  {
    for (int j = 0; j < NUM_SEAVEYS * 2; j++) 
    {
      int ant = j /2; 
      int pol = j %2; 

      TH2F * h = isRF ? avgs[ant][pol] : avgs_minbias[ant][pol]; 
      TH1I * hnorm = isRF ? norms : norms_minbias; 


      for (int ybin = 1; ybin <= h->GetNbinsY(); ybin++)
      {
        if (hnorm->GetBinContent(ybin) < min_norm) 
        {
          continue; 
        }

        for (int k = 0; k<h->GetNbinsX() ; k++)
        {
            h->SetBinContent(k,ybin, h->GetBinContent(k,ybin) / hnorm->GetBinContent(ybin)); 
        }

        if (!isRF) 
        {
            rms[ant][pol]->SetBinContent(ybin, rms[ant][pol]->GetBinContent(ybin) / hnorm->GetBinContent(ybin)); 
        }

      }

      //now go through and fill in the rows with no content with an average of ones that do 
      for (int ybin = 1; ybin <= h->GetNbinsY(); ybin++)
      {

        if (hnorm->GetBinContent(ybin) < min_norm)
        {
          int last_bin = ybin-1; 
          while(hnorm->GetBinContent(last_bin) < min_norm && last_bin > 0) last_bin--; 

          int next_bin = ybin+1; 
          while (!hnorm->GetBinContent(next_bin) && next_bin <= h->GetNbinsY()) next_bin++; 

          //lost cause, no good bins
          if (last_bin == 0 && next_bin > h->GetNbinsY()) break; 

          if (last_bin == 0) 
          {
            for (int k = 0; k<h->GetNbinsX() ; k++)
            {
                h->SetBinContent(k,ybin, h->GetBinContent(k,next_bin)); 
                if (!isRF) rms[ant][pol]->SetBinContent(ybin, rms[ant][pol]->GetBinContent(next_bin)); 
            }
          }
          else if (next_bin > h->GetNbinsY()) 
          {
            for (int k = 0; k<h->GetNbinsX() ; k++)
            {
                h->SetBinContent(k,ybin, h->GetBinContent(k,last_bin)); 
                if (!isRF) rms[ant][pol]->SetBinContent(ybin, rms[ant][pol]->GetBinContent(last_bin)); 
            }
          }
          else
          {
            double last_frac = double(ybin - last_bin) / (next_bin-last_bin); 
            
            for (int k = 0; k<h->GetNbinsX() ; k++)
            {
                h->SetBinContent(k,ybin, last_frac * h->GetBinContent(k,last_bin)  + (1-last_frac) * h->GetBinContent(k,next_bin)); 
                if (!isRF) rms[ant][pol]->SetBinContent(ybin, last_frac * rms[ant][pol]->GetBinContent(last_bin) + (1-last_frac) * rms[ant][pol]->GetBinContent(next_bin)); 
            }
          }
        }
      }
    }

  }

  return 0; 
}



void UCorrelator::TimeDependentAverage::saveToDir(const char * dir) 
{

  gSystem->mkdir(dir,true); 

  TString fstr; 
  fstr.Form("%s/%d_%d.root", dir, run,nsecs); 
  TFile f(fstr.Data(),"RECREATE"); 
  f.cd(); 

  f.mkdir("specavg"); 
  f.cd("specavg"); 
  for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
  {
     avgs[ant][0]->Write(TString::Format("h%d",ant)); 
     avgs[ant][1]->Write(TString::Format("v%d",ant)); 
  }

  f.cd(); 

  f.mkdir("specavg_minbias"); 
  f.cd("specavg_minbias"); 
  for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
  {
     avgs_minbias[ant][0]->Write(TString::Format("h%d",ant)); 
     avgs_minbias[ant][1]->Write(TString::Format("v%d",ant)); 
  }
  f.cd();

  f.mkdir("rms"); 
  f.cd("rms"); 
  for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
  {
     rms[ant][0]->Write(TString::Format("h%d",ant)); 
     rms[ant][1]->Write(TString::Format("v%d",ant)); 
  }
  f.cd();

  nblasts->Write("nblasts"); 
  norms->Write("norms"); 
  norms_minbias->Write("norms_minbias"); 


}

UCorrelator::TimeDependentAverage::TimeDependentAverage(int run, int nsecs, const char * persistdir,
      double max_bottom_top_ratio, int min_norm, double max_power) 
  : nsecs(nsecs) , run(run)
{
  bool foundit = false;

  const char * check_dir = persistdir ? persistdir 
                                      :  getenv("UCORRELATOR_TIMEAVG_DIR") ? : 0; 


  if (check_dir)
  {

    TFile f(TString::Format("%s/%d_%d.root", check_dir, run,nsecs)); 
    if (f.IsOpen())
    {
      f.cd(); 


      TH1I * found= (TH1I*) gDirectory->Get("norms"); 
      norms = new TH1I(*found); 
      norms->SetDirectory(0); 

      found= (TH1I*) gDirectory->Get("norms_minbias"); 
      norms_minbias = new TH1I(*found); 
      norms_minbias->SetDirectory(0); 

      found= (TH1I*) gDirectory->Get("nblasts"); 
      nblasts = new TH1I(*found); 
      nblasts->SetDirectory(0); 

      f.cd("specavg"); 
      for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
      {
        TH2F * found_hpol = (TH2F*) gDirectory->Get(TString::Format("h%d",ant)); 
        avgs[ant][0] = new TH2F(*found_hpol);  
        avgs[ant][0]->SetDirectory(0); 

        TH2F * found_vpol = (TH2F*) gDirectory->Get(TString::Format("v%d",ant)); 
        avgs[ant][1] = new TH2F(*found_vpol);  
        avgs[ant][1]->SetDirectory(0); 
      }

      f.cd("specavg_minbias"); 

      for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
      {

        TH2F * found_hpol = (TH2F*) gDirectory->Get(TString::Format("h%d",ant)); 
        avgs_minbias[ant][0] = new TH2F(*found_hpol);  
        avgs_minbias[ant][0]->SetDirectory(0); 

        TH2F * found_vpol = (TH2F*) gDirectory->Get(TString::Format("v%d",ant)); 
        avgs_minbias[ant][1] = new TH2F(*found_vpol);  
        avgs_minbias[ant][1]->SetDirectory(0); 
      }

      f.cd("rms"); 

      for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
      {

        TH1D * found_hpol = (TH1D*) gDirectory->Get(TString::Format("h%d",ant)); 
        rms[ant][0] = new TH1D(*found_hpol);  
        rms[ant][0]->SetDirectory(0); 

        TH1D * found_vpol = (TH1D*) gDirectory->Get(TString::Format("v%d",ant)); 
        rms[ant][1] = new TH1D(*found_vpol);  
        rms[ant][1]->SetDirectory(0); 
      }

      foundit = true;
    }
  }

  if (!foundit) 
  {
    printf("Didn't find spectrum_averages... creating!\n"); 
    if (!check_dir) 
    {
      printf("Define UCORRELATOR_TIMEAVG_DIR to persist somewhere.\n"); 
    }
    computeAverage(  max_bottom_top_ratio, min_norm, max_power); 
    if (check_dir) saveToDir(check_dir); 
  }

  memset(peakiness,0,sizeof(peakiness)); 
  memset(peakiness_minbias,0,sizeof(peakiness_minbias)); 
}
            


TH1* UCorrelator::TimeDependentAverage::getSpectrumAverage(AnitaPol::AnitaPol_t pol, int ant, double t, bool db, bool minbias)  const
{

  const TH2 * h = minbias ? avgs_minbias[ant][pol] : avgs[ant][pol];
  //figure out which bin we are in 
   int bin =  h->GetYaxis()->FindFixBin(t); 

   if (bin == 0 || bin == h->GetNbinsY()+1) return 0 ; 

   TString title; 
   title.Form("Ant%d %c-Pol Average at t=%g", ant, pol == AnitaPol::kHorizontal ? 'H' : 'V', t); 
   TString name;
   name.Form("specavg_%d_%d_%g",ant,pol,t); 
   TH1D * answer =  new TH1D(name,title, 131,0,1.31); //todo don't hardcode 

   for (int i = 1; i <= answer->GetNbinsX(); i++)
   {
     answer->SetBinContent(i, h->GetBinContent(i,bin)); 
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

TH1 *UCorrelator::TimeDependentAverage::getSpectrumPercentile(AnitaPol::AnitaPol_t pol, int ant, double pct , bool db, bool minbias ) const
{

  TH1 * answer = UCorrelator::image::getPctileProjection( minbias ? avgs_minbias[ant][pol] : avgs[ant][pol], 1, pct, true, minbias ? norms_minbias : norms); 

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


double UCorrelator::TimeDependentAverage::getStartTime() const 
{

  return avgs[0][0] ? avgs[0][0]->GetYaxis()->GetXmin() : -1; 
}

double UCorrelator::TimeDependentAverage::getEndTime() const 
{
  return avgs[0][0] ? avgs[0][0]->GetYaxis()->GetXmax() : -1;   
}


static const UCorrelator::TimeDependentAverage* defaultThermalAvg = 0; 
const UCorrelator::TimeDependentAverage* UCorrelator::TimeDependentAverage::defaultThermal() 
{
  if (defaultThermalAvg) return defaultThermalAvg; 

  //TODO make this a4 compatible as well 
   if (AnitaVersion::get() == 4) 
   {
     fprintf(stderr,"warning: using default terminated thermal spectrum for A3 for peakiness\n"); 
   }

   TString dir; 
   dir.Form("%s/share/UCorrelator/terminated_noise/", getenv("ANITA_UTIL_INSTALL_DIR")); 

   defaultThermalAvg = new TimeDependentAverage(11382,60, dir.Data()); 
   return defaultThermalAvg; 
}



void UCorrelator::TimeDependentAverage::computePeakiness(const TimeDependentAverage * thermalSpec, double fractionForNormalization) const
{

  if (!thermalSpec) thermalSpec = defaultThermal(); 


//#ifdef UCORRELATOR_OPENMP
//#pragma omp parallel for
//#endif 
  for (int ant  = 0; ant < NUM_SEAVEYS; ant++)
  {
//    printf("%d\n",ant);
    for (int ipol = 0; ipol < 2; ipol++) 
    {

      for (int minbias = 0; minbias < 2; minbias++)
      {
//      printf("%d %d\n",ant,ipol);

        //get median spectrum 
        TH1 * median = getSpectrumPercentile(AnitaPol::AnitaPol_t(ipol), ant, 0.5,false,minbias); 
        TH1 * thermal = thermalSpec->getSpectrumPercentile(AnitaPol::AnitaPol_t(ipol),ant,0.5,false,minbias);

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
        name.Form("peakiness_%d_%d_%s\n", ant,ipol, minbias ? "minbias" : "rf"); 
        TString title; 
        title.Form("%s peakiness ant=%d pol=%d\n", minbias ? "minbias" : "RF" , ant,ipol); 

        TH2 * avg = minbias ? avgs_minbias[ant][ipol] : avgs[ant][ipol]; 
        TH2D * peaky = new TH2D(name,title,
                                         avg->GetNbinsX(), avg->GetXaxis()->GetXmin(), avg->GetXaxis()->GetXmax(), 
                                         avg->GetNbinsY(), avg->GetYaxis()->GetXmin(), avg->GetYaxis()->GetXmax());  


        peaky->SetDirectory(0); 

        for (int ii = 1; ii < peaky->GetNbinsX(); ii++)
        {
          for (int jj = 1; jj < peaky->GetNbinsY(); jj++)
          {
            peaky->SetBinContent(ii,jj, avg->GetBinContent(ii,jj)  / ( ratio * thermal->GetBinContent(ii))); 
            if (std::isnan(peaky->GetBinContent(ii,jj)))
                peaky->SetBinContent(ii,jj,0); 
          }
        }
         (minbias ? peakiness_minbias[ant][ipol] : peakiness[ant][ipol]) = peaky; 
        delete median; 
        delete thermal; 
      }

    }
  }

}

UCorrelator:: TimeDependentAverageLoader::TimeDependentAverageLoader(const char * the_dir, int secs) 
  : dir(the_dir), nsecs(secs)
{
  tavg = 0; 
}

static TMutex m; 

const UCorrelator::TimeDependentAverage* UCorrelator::TimeDependentAverageLoader::avg(double t) const
{

  m.Lock(); 
//  printf("%g\n",t); 

  if (tavg && t >= tavg->getStartTime()-5 && t <= tavg->getEndTime()+5 ) 
  {
    m.UnLock(); 
    return tavg; 
  }

  int run = AnitaDataset::getRunAtTime(t); 
  printf("loading average from run %d\n",run); 
  
  if (tavg) delete tavg; 
  tavg = new TimeDependentAverage(run,nsecs, dir); 
  tavg->computePeakiness(); 
  m.UnLock(); 
  
  return tavg; 

}

double UCorrelator::TimeDependentAverage::getBlastFraction(double t) const
{
  //Estimate by estimating nblasts and n, taking their ratio
  int bin = nblasts->GetXaxis()->FindFixBin(t); 

  //the other closest bin
  int other_bin = t < nblasts->GetXaxis()->GetBinCenter(bin) ? bin-1 : bin+1; 


  //the fraction of blasts in the closest bin
  double blast_frac_bin = nblasts->GetBinContent(bin) / (nblasts->GetBinContent(bin) + norms->GetBinContent(bin)); 

  //f is the weight of the main bin. 
  double f = norms->GetBinContent(other_bin) == 0 ? 1 :  1-fabs(t-nblasts->GetXaxis()->GetBinCenter(bin)) / nsecs; 
  if (f < 1 )
  {
    double blast_frac_other_bin = nblasts->GetBinContent(other_bin) / (nblasts->GetBinContent(other_bin) + norms->GetBinContent(other_bin)); 
    return f * blast_frac_bin + (1-f) * blast_frac_other_bin; 
  }

  return blast_frac_bin; 
}



static UCorrelator::TimeDependentAverageLoader * getLoader(int nsecs) 
{
  static __thread UCorrelator::TimeDependentAverageLoader * ldr = 0; 
  if (!ldr || ldr->getNsecs() != nsecs) 
  {
    if (ldr) delete ldr; 
    ldr = new UCorrelator::TimeDependentAverageLoader(0,nsecs); 
  }

  return ldr; 
}

double UCorrelator::TimeDependentAverage::getRMS(AnitaPol::AnitaPol_t pol, int ant, double t) const
{
  return ((TH1*) getRMS(pol,ant))->Interpolate(t); 
}



double UCorrelator::TimeDependentAverageLoader::getRMS(double t, AnitaPol::AnitaPol_t pol, int ant, int nsecs) 
{
  return getLoader(nsecs)->avg(t)->getRMS(pol,ant,t); 
}

double UCorrelator::TimeDependentAverageLoader::getPayloadBlastFraction(double t,  int nsecs) 
{
  return getLoader(nsecs)->avg(t)->getBlastFraction(t); 
}


const TH2D * UCorrelator::TimeDependentAverage::getPeakiness(AnitaPol::AnitaPol_t pol, int ant, bool minbias) const
{

  m.Lock();
  if (!peakiness[0][0])
  {
    computePeakiness(); 
  }
  m.UnLock(); 

  return minbias ? peakiness_minbias[ant][pol] : peakiness[ant][pol]; 
}
