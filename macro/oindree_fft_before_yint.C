//int N = 262144; 
//int N = 128; 

void oindree_fft_before_yint(int whichrun, int whichanita, const int N = 262144);

void oindree_fft_before_yint(int whichrun, int whichanita, const int N)

{
  cout << "Beginning... " << endl; 

  AnitaDataset d(whichrun); 
  AnitaVersion::set(whichanita); 

  cout << "Number of entries in this run is " << d.N() << endl; 

  d.getEntry(0);
  int start_time = d.header()->realTime; 
  cout << "realTime is " << start_time << endl;
  
  d.getEntry(d.N()-1);
  int end_time = d.header()->realTime;
  cout << "realTime is " << end_time << endl;

  TH1D hrealTime("hrealTime",Form("ANITA-%i Data Run %i;Time (Month-Day);Number of Events",whichanita,whichrun),N,start_time,end_time);

  ////to check for events just before y-int cut
  TChain beforeyint("analysisOutput0");
  beforeyint.Add("/data/anita/oindree/anita3analysisoutput/fromJacob/analysisOutput_2_175_439_CpolInfo.root");
  cout << beforeyint.GetEntries() << endl; 

  TCanvas c("c","c",1000,800); 
  hrealTime.GetXaxis()->SetTimeDisplay(1); 
  hrealTime.GetXaxis()->SetTimeFormat("%m-%d");
  hrealTime.GetXaxis()->SetNdivisions(-307);  
  beforeyint.Draw("realTime >> hrealTime",Form("runNumber == %i",whichrun),"");
  //beforeyint.Draw("realTime >> hrealTime","","");
  c.SetLogy(); 
  c.SaveAs(Form("realTimeRun%iAnita%iBeforeYint.png",whichrun,whichanita));
  c.SaveAs(Form("realTimeRun%iAnita%iBeforeYint.root",whichrun,whichanita));

  ///////////////////////////////////////////////////// FFT ///////////////////////////////////////////////////////////////////


  cout << "Beginning FFT business ..." << endl; 

  double y[N];
  double realtimes[N]; 
  double newy[N/2]; 
  double freqs[N/2];  

  for (int ibin = 0; ibin < N; ibin++)
    {
      y[ibin] = hrealTime.GetBinContent(ibin); 
      //cout << hrealTimeJan8.GetBinContent(ibin) << endl;
      //cout << "hrealTime getbincenter " << setprecision(10) << hrealTime.GetBinCenter(ibin) << endl;
      realtimes[ibin] = hrealTime.GetBinCenter(ibin) - hrealTime.GetBinCenter(0); 
      //cout << "realtimes array " << realtimes[ibin] << endl; 
      //cout << ibin << endl;
    } 

  TGraph gr(N,realtimes,y); 
  TCanvas yy("yy","yy",1000,800);
  gr.SetTitle("the thing I will fft"); 
  gr.GetXaxis()->SetTitle("Time (second)");
  gr.GetYaxis()->SetTitle("Number of Events"); 
  //gr_original.GetYaxis()->SetRangeUser(0,100); 
  gr.Draw("apl");
  yy.SetLogy(); 
  yy.SaveAs(Form("realTimeRun%iAnita%iOriginalBeforeYint.png",whichrun,whichanita)); 


  FFTWComplex *fft = FFTtools::doFFT(N,y);
  
  double deltaT = (realtimes[N-1] - realtimes[0])/N; //realtimes is in seconds 
  double deltaF = 1/(N * deltaT); //will be in Hz
  cout << "N is " << N << endl; 
  cout << "deltaT is " << deltaT << endl; 
  cout << "deltaF is " << setprecision(3) << deltaF << endl; 

  double tempF = deltaF; 
 
  for (int j = 0; j < (N/2)+1; j++) 
    {
      newy[j] = FFTtools::getAbs(fft[j]); 
      freqs[j] = tempF;
      tempF = tempF + deltaF; 
    }
  
  delete [] fft;

  TLine line(1,1e-1,1,1e5); 
  line.SetLineColor(2); 

  TGraph gr_fft(N/2,freqs,newy); 
  TCanvas f("f","f",1000,800); 
  gr_fft.SetTitle(Form("fft:N=%i,dF=%fuHz",N,deltaF*1e6)); 
  gr_fft.GetXaxis()->SetTitle("Frequency (Hz)");
  gr_fft.GetYaxis()->SetTitle("FFT of Number of Events"); 
  //gr_fft.GetXaxis()->SetRangeUser(deltaF,5*deltaF); 
  gr_fft.Draw("apl");
  line.Draw("same"); 
  f.SetLogy();
  //gr_fft.GetXaxis()->SetRangeUser(deltaF,(N/2)*deltaF); 
  f.SetLogx();  
  f.SaveAs(Form("realTimeRun%iAnita%iFFTBeforeYint.png",whichrun,whichanita));
  f.SaveAs(Form("realTimeRun%iAnita%iFFTBeforeYint.root",whichrun,whichanita)); 

  cout << gr_fft.Eval(1.) << endl; 
  double amplitude = gr_fft.Eval(1.); 
  
  double ywhosexis1 = sqrt( pow(amplitude,2) - 1 ); 

  cout << ywhosexis1 << endl; 
  double phase = TMath::ATan2(ywhosexis1,1.); 

  cout << "phase is " << phase << endl; 
  













}
