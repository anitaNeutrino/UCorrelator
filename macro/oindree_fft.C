//int N = 262144; 
//int N = 128; 

void oindree_fft(int whichrun, int whichanita, const int N = 262144);

void oindree_fft(int whichrun, int whichanita, const int N)

{

  cout << "Beginning... " << endl; 

  AnitaDataset d(whichrun); 
  AnitaVersion::set(whichanita); 

  cout << "Number of entries in this run is " << d.N() << endl; 

  d.getEntry(0);
  double start_time = d.header()->triggerTime + (d.header()->triggerTimeNs)/1e9; 
  cout << "start time is " << start_time << endl;
  cout << "trigger time ns of the first entry " << setprecision(10) << d.header()->triggerTimeNs << endl; 
  
  d.getEntry(d.N()-1);
  double end_time = (d.header()->triggerTime + (d.header()->triggerTimeNs)/1e9) - start_time;
  cout << "end time is " << end_time << endl;

  UsefulAdu5Pat pat(d.gps());  

  TH1D hrealTime("hrealTime",Form("ANITA-%i Data Run %i;Time (Month-Day);Number of Events",whichanita,whichrun),N,0,end_time);

  double time = 0.; 

  for (int ientry = 0; ientry < d.N(); ientry++) {
    d.getEntry(ientry);
    //cout << d.header()->trigTypeAsString() << endl; 
    if ((strcmp(d.header()->trigTypeAsString(),"RF") == 0) && (UCorrelator::isLDB(d.header()) == 0) && (UCorrelator::isWAISHPol(&pat,d.header()) == 0) && (UCorrelator::isWAISVPol(&pat,d.header()) == 0)) 
    {   /// if string compare of the thing and "RF" is equal to 0 then they are the same string
      time = (d.header()->triggerTime - start_time) + (d.header()->triggerTimeNs)/1e9; 
      //hrealTime.Fill(d.header()->realTime);
      hrealTime.Fill(time); 
      //cout << "wais hpol " << UCorrelator::isWAISHPol(&pat,d.header()) << endl;
      //cout << "wais vpol " << UCorrelator::isWAISVPol(&pat,d.header()) << endl; 
    }
    else continue; 
  }

  TCanvas c("c","c",1000,800); 
  //hrealTime.GetXaxis()->SetTimeDisplay(1); 
  //hrealTime.GetXaxis()->SetTimeFormat("%m-%d");
  //hrealTime.GetXaxis()->SetNdivisions(-307);  
  hrealTime.Draw();  // commented out to check for events right before y-int cut
  c.SetLogy(); 
  c.SaveAs(Form("realTimeRun%iAnita%i.png",whichrun,whichanita));
  c.SaveAs(Form("realTimeRun%iAnita%i.root",whichrun,whichanita));

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
      realtimes[ibin] = hrealTime.GetBinCenter(ibin); 
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
  yy.SaveAs(Form("realTimeRun%iAnita%iOriginal.png",whichrun,whichanita)); 


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

  TLine line(1,1e-2,1,1e5); 
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
  f.SaveAs(Form("realTimeRun%iAnita%iFFT.png",whichrun,whichanita));
  f.SaveAs(Form("realTimeRun%iAnita%iFFT.root",whichrun,whichanita)); 

  cout << gr_fft.Eval(1.) << endl; 
  double amplitude = gr_fft.Eval(1.); 
  
  double ywhosexis1 = sqrt( pow(amplitude,2) - 1 ); 

  cout << ywhosexis1 << endl; 
  double phase = TMath::ATan2(ywhosexis1,1.); 

  cout << "phase is " << phase << endl; 

  













}
