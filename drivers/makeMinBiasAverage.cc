#include "FFTtools.h"
#include "TH2.h" 
#include "FilterStrategy.h" 
#include "Analyzer.h" 
#include "UCFilters.h" 
#include "TFile.h" 
#include "AnitaDataset.h" 
#include "UsefulAdu5Pat.h" 
#include "FilteredAnitaEvent.h" 
#include <stdio.h> 
#include <stdlib.h> 
#include "RawAnitaHeader.h" 
#include "TTree.h" 

void makeMinBiasAverage(int run, int max = 0, int partial_time = 60, int filter = 0, const char * outdir = "mbavg", int anita=4) 
{


  FFTtools::loadWisdom("wisdom.dat"); 
  FilterStrategy strat; 
  if (filter) 
    UCorrelator::fillStrategyWithKey(&strat,"sinsub_05_3_ad_2"); 


  AnitaDataset d(run,false, WaveCalType::kDefault, (AnitaDataset::DataDirectory) anita, AnitaDataset::kDefault );

  UCorrelator::Correlator corr_h(900, 0,360,450,-45,45,true); 
  UCorrelator::Correlator corr_v(900, 0,360,450,-45,45,true); 

  TString mkdir = Form("mkdir -p %s/a%d", outdir,anita); 
  system(mkdir.Data()); 
  TFile of(Form("%s/a%d/avg%d%s.root",outdir, anita,run,filter?"_filtered":""), "RECREATE");


  TH2D * equatorial_h = new TH2D("equatorial_h","HPol (Equatorial) ; RA (hr); dec(deg)", 1800, 0,24,900, -90,90); 
  TH2D * equatorial_v = new TH2D("equatorial_v","VPol (Equatorial) ; RA (hr); dec(deg)", 1800, 0,24,900, -90,90); 
  TH2D * equatorial_norm = new TH2D("equatorial_norm","Normalization (Equatorial) ; RA (hr); dec(deg)", 1800, 0,24,900, -90,90); 

  TH2D * sun_h = new TH2D("sun_h","HPol (Sun-Centered) ; #Delta phi (deg); #Delta theta (deg)", 200, -10,10,200,-10,10); 
  TH2D * sun_v = new TH2D("sun_v","VPol (Sun-Centered) ; #Delta phi (deg); #Delta theta (deg)", 200, -10,10,200,-10,10); 
  TH2D * sun_norm = new TH2D("sun_norm","Normalization (Sun-Centered) ; #Delta phi (deg); #Delta theta (deg)", 200, -10,10,200,-10,10); 

  TH2D * horiz_h = new TH2D("horiz_h","HPol (Horizontal Coords) ; bearing (deg); elevation above horiz (deg)", 1800, 0,360,900,-45,45); 
  TH2D * horiz_v = new TH2D("horiz_v","VPol (Horizontal Coords) ; bearing (deg); elevation above horiz (deg)", 1800, 0,360,900,-45,45); 
  TH2D * horiz_norm = new TH2D("horiz_norm","Normalization (Horizontal Coords) ; bearing (deg); elevation above horiz (deg)", 1800, 0,360,900,-45,45); 




  TTree * t = new TTree("time_dependent","Time Dependent maps"); 

  int start_time = 0; 
  int end_time = 0;

  double sun_ra;
  double sun_dec;
  double sun_max_el_h;
  double sun_max_br_h;
  double sun_max_el_v;
  double sun_max_br_v;
  double sun_el; 
  double sun_br; 
  double sun_anita_phi; 
  t->Branch("start_time",&start_time);
  t->Branch("end_time",&end_time);
  t->Branch("sun_ra",&sun_ra);
  t->Branch("sun_dec",&sun_dec);
  t->Branch("sun_el",&sun_el);
  t->Branch("sun_br",&sun_br);
  t->Branch("sun_max_el_h",&sun_max_el_h);
  t->Branch("sun_max_br_h",&sun_max_br_h);
  t->Branch("sun_max_el_v",&sun_max_el_v);
  t->Branch("sun_max_br_v",&sun_max_br_v);
  t->Branch("sun_anita_phi",&sun_anita_phi);


  TH2D * equatorial_h_partial = (TH2D*) equatorial_h->Clone("equatorial_h_partial"); 
  TH2D * equatorial_v_partial = (TH2D*) equatorial_v->Clone("equatorial_v_partial"); 
  TH2D * equatorial_norm_partial = (TH2D*) equatorial_norm->Clone("equatorial_norm_partial"); 

  TH2D * sun_h_partial = (TH2D*) sun_h->Clone("sun_h_partial"); 
  TH2D * sun_v_partial = (TH2D*) sun_v->Clone("sun_v_partial"); 
  TH2D * sun_norm_partial = (TH2D*) sun_norm->Clone("sun_norm_partial"); 

  TH2D * horiz_h_partial = (TH2D*) horiz_h->Clone("horiz_h_partial"); 
  TH2D * horiz_v_partial = (TH2D*) horiz_v->Clone("horiz_v_partial"); 
  TH2D * horiz_norm_partial = (TH2D*) horiz_norm->Clone("horiz_norm_partial"); 


  t->Branch("equatorial_h",&equatorial_h_partial); 
  t->Branch("equatorial_v",&equatorial_v_partial); 
  t->Branch("equatorial_norm",&equatorial_norm_partial); 
  t->Branch("sun_h",&sun_h_partial); 
  t->Branch("sun_v",&sun_v_partial); 
  t->Branch("sun_norm",&sun_norm_partial); 
  t->Branch("horiz_h",&horiz_h_partial); 
  t->Branch("horiz_v",&horiz_v_partial); 
  t->Branch("horiz_norm",&horiz_norm_partial); 




  int ndone = 0;
 
  for (int i =0; i < d.N(); i++) 
  {
    d.getEntry(i); 
    if (d.header()->getTriggerBitRF()) continue; 
    printf("%d\n", i); 

    if (start_time ==0) start_time = d.header()->triggerTime; 

    FilteredAnitaEvent ev(d.useful(), &strat, d.gps(), d.header()); 

    corr_v.compute(&ev, AnitaPol::kVertical); 
    corr_h.compute(&ev, AnitaPol::kHorizontal); 

    UsefulAdu5Pat pat(d.gps()); 


    TH2 * h = (TH2*) corr_h.getHist();
    TH2 * v = (TH2*) corr_v.getHist();



    int nx = equatorial_h->GetNbinsX(); 
    int ny = equatorial_h->GetNbinsY(); 
    double xmin = equatorial_h->GetXaxis()->GetXmin();
    double ymin = equatorial_h->GetYaxis()->GetXmin();
    double dx = equatorial_h->GetXaxis()->GetBinWidth(1);
    double dy = equatorial_h->GetYaxis()->GetBinWidth(1);


   //equatorial 

#pragma omp parallel for
    for (int ra = 1; ra <= nx; ra++)
    {
      double x=  (ra-0.5)*dx+xmin;

      double v_phi[ny];
      double v_theta[ny];
      double v_ra[ny];
      double v_dec[ny];
      int v_bin[ny];
      int ii = 0; 

      for (int dec = 1; dec <= ny; dec++)
      {
        double y= (dec-0.5)*dy + ymin; 
        v_ra[ii] = x; 
        v_dec[ii] = y; 
        int ibin = ra + dec * (nx+2);
        v_bin[ii] = ibin; 
        ii++;
      }

      pat.fromRADec(ny, &v_ra[0],&v_dec[0], &v_phi[0], &v_theta[0]);

      for (ii = 0; ii < ny; ii++)
      {
          if (v_theta[ii] > -45 && v_theta[ii] < 5) 
          {
            equatorial_h->GetArray()[v_bin[ii]] += h->Interpolate(v_phi[ii], -v_theta[ii]);
            equatorial_v->GetArray()[v_bin[ii]] += v->Interpolate(v_phi[ii], -v_theta[ii]);
            equatorial_norm->GetArray()[v_bin[ii]]++; 
            equatorial_h_partial->GetArray()[v_bin[ii]] += h->Interpolate(v_phi[ii], -v_theta[ii]);
            equatorial_v_partial->GetArray()[v_bin[ii]] += v->Interpolate(v_phi[ii], -v_theta[ii]);
            equatorial_norm_partial->GetArray()[v_bin[ii]]++; 
          }

      }
    }

   

    //horizontal
    nx = horiz_h->GetNbinsX(); 
    ny = horiz_h->GetNbinsY(); 
    xmin = horiz_h->GetXaxis()->GetXmin();
    ymin = horiz_h->GetYaxis()->GetXmin();
    dx = horiz_h->GetXaxis()->GetBinWidth(1);
    dy = horiz_h->GetYaxis()->GetBinWidth(1);


#pragma omp parallel for
    for (int br = 1; br <= nx; br++)
    {
      double phi,theta; 
      double x=  (br-0.5) * dx + xmin;
      for (int el = 1; el <= ny; el++)
      {
        double y= (el-0.5)*dy +ymin;
        int ibin = br + el * (nx+2);
        theta = -y;
        phi = FFTtools::wrap(pat.heading-x,360);

        horiz_h->GetArray()[ibin] += h->Interpolate(phi,-theta); 
        horiz_v->GetArray()[ibin] += v->Interpolate(phi,-theta); 
        horiz_norm->GetArray()[ibin]++;

        horiz_h_partial->GetArray()[ibin]+= h->Interpolate(phi,-theta); 
        horiz_v_partial->GetArray()[ibin]+= v->Interpolate(phi,-theta); 
        horiz_norm_partial->GetArray()[ibin]++;


      }
    }


    //sun 
    nx = sun_h->GetNbinsX(); 
    ny = sun_h->GetNbinsY(); 
    xmin = sun_h->GetXaxis()->GetXmin();
    ymin = sun_h->GetYaxis()->GetXmin();
    dx = sun_h->GetXaxis()->GetBinWidth(1);
    dy = sun_h->GetYaxis()->GetBinWidth(1);


    double sun_phi, sun_theta;
    pat.getSunPosition(sun_phi,sun_theta);
#pragma omp parallel for
    for (int dphi = 1; dphi <= nx; dphi++)
    {
      double phi,theta; 
      double x = xmin + (dphi-0.5) * dx; 
      for (int dtheta =1; dtheta <= ny; dtheta++)
      {
        double y = ymin + (dtheta-0.5)*dy;

        phi = FFTtools::wrap(sun_phi+x,360);
        theta = sun_theta-y;
        

        if ( theta > -45 && theta < 5)//avoid the horizon, for now
        {
          int ibin = dphi + dtheta * (nx+2);

          sun_h->GetArray()[ibin]+= h->Interpolate(phi,-theta); 
          sun_v->GetArray()[ibin]+= v->Interpolate(phi,-theta); 
          sun_norm->GetArray()[ibin]++; 

          sun_h_partial->GetArray()[ibin]+= h->Interpolate(phi,-theta); 
          sun_v_partial->GetArray()[ibin]+= v->Interpolate(phi,-theta); 
          sun_norm_partial->GetArray()[ibin]++; 


        }
      }
    }

    int last_iteration = (max > 0 && ++ndone > max) || i == d.N()-1; 

    if (((int) d.header()->triggerTime) - start_time > partial_time || last_iteration) 
    {


      of.cd(); 

      equatorial_h_partial->Divide(equatorial_norm_partial);
      equatorial_v_partial->Divide(equatorial_norm_partial);
      horiz_h_partial->Divide(horiz_norm_partial);
      horiz_v_partial->Divide(horiz_norm_partial);
      sun_h_partial->Divide(sun_norm_partial);
      sun_v_partial->Divide(sun_norm_partial);

      equatorial_h_partial->SetEntries(equatorial_norm_partial->Integral());
      equatorial_v_partial->SetEntries(equatorial_norm_partial->Integral());
      equatorial_norm_partial->SetEntries(equatorial_norm_partial->Integral());
      horiz_h_partial->SetEntries(horiz_norm_partial->Integral());
      horiz_v_partial->SetEntries(horiz_norm_partial->Integral());
      horiz_norm_partial->SetEntries(horiz_norm_partial->Integral());
      sun_h_partial->SetEntries(sun_norm_partial->Integral());
      sun_v_partial->SetEntries(sun_norm_partial->Integral());
      sun_norm_partial->SetEntries(sun_norm_partial->Integral());

      end_time = d.header()->triggerTime; 

      sun_br = FFTtools::wrap(pat.heading-sun_phi,360);
      sun_el = -sun_theta;
      pat.astronomicalCoordinates(sun_phi,sun_theta,&sun_ra, &sun_dec); 

      int max_bin_h = sun_h_partial->GetMaximumBin(); 
      int max_bin_v = sun_v_partial->GetMaximumBin(); 
      int bx,by,bz; 
      sun_h_partial->GetBinXYZ(max_bin_h, bx,by,bz); 
      sun_max_el_h =sun_h_partial->GetYaxis()->GetBinCenter(by); 
      sun_max_br_h =sun_h_partial->GetXaxis()->GetBinCenter(bx); 
      sun_v_partial->GetBinXYZ(max_bin_v, bx,by,bz); 
      sun_max_el_v =sun_v_partial->GetYaxis()->GetBinCenter(by); 
      sun_max_br_v =sun_v_partial->GetXaxis()->GetBinCenter(bx); 

      sun_anita_phi = sun_phi;


      t->Fill();

      start_time = end_time; 
      equatorial_h_partial->Reset();
      equatorial_v_partial->Reset();
      equatorial_norm_partial->Reset();
      horiz_h_partial->Reset();
      horiz_v_partial->Reset();
      horiz_norm_partial->Reset();
      sun_h_partial->Reset();
      sun_v_partial->Reset();
      sun_norm_partial->Reset();


    }


    if (last_iteration) break; 
  }

  equatorial_h->Divide(equatorial_norm);
  equatorial_v->Divide(equatorial_norm);
  horiz_h->Divide(horiz_norm);
  horiz_v->Divide(horiz_norm);
  sun_h->Divide(sun_norm);
  sun_v->Divide(sun_norm);

  equatorial_h->SetEntries(equatorial_norm->Integral());
  equatorial_v->SetEntries(equatorial_norm->Integral());
  equatorial_norm->SetEntries(equatorial_norm->Integral());
  horiz_h->SetEntries(horiz_norm->Integral());
  horiz_v->SetEntries(horiz_norm->Integral());
  horiz_norm->SetEntries(horiz_norm->Integral());
  sun_h->SetEntries(sun_norm->Integral());
  sun_v->SetEntries(sun_norm->Integral());
  sun_norm->SetEntries(sun_norm->Integral());



  equatorial_h->Write();
  equatorial_v->Write();
  equatorial_norm->Write();
  horiz_h->Write();
  horiz_v->Write();
  horiz_norm->Write();
  sun_h->Write();
  sun_v->Write();
  sun_norm->Write();

  t->Write(); 


  FFTtools::saveWisdom("wisdom.dat"); 

}

int main(int nargs, char ** args) 
{
  int run = atoi(args[1]); 
  int max = nargs > 2 ?atoi(args[2]) : 0; 
  int partial_time = nargs > 3 ?atoi(args[3]) : 60; 
  int filter = nargs > 4 ?atoi(args[4]) : 0; 
  int anita = nargs > 5 ? atoi(args[5]) : 4; 
  const char * outdir = nargs > 6 ? args[6] :"mbavg";
  makeMinBiasAverage(run,max,partial_time,filter, outdir,anita); 

}
