void plotSatellites(const char * file, const char * tree="decimated_no_filter") 
{ 


  TChain c(tree); 
  c.Add(file); 



  TCanvas * cv = new TCanvas; 

  cv->Divide(4,2); 


  cv->cd(1); 
  c.Draw("coherent[0][0].I", "peak[0][0].theta <-10 && peak[0][0].snr > 3"); 
  cv->cd(2); 
  c.Draw("coherent[0][0].Q / coherent[0][0].I", "peak[0][0].theta <-10 && peak[0][0].snr > 3 && coherent[0][0].I > 0"); 
  cv->cd(3); 
  c.Draw("coherent[0][0].U / coherent[0][0].I", "peak[0][0].theta <-10 && peak[0][0].snr > 3"); 
  cv->cd(4); 
  c.Draw("coherent[0][0].V / coherent[0][0].I", "peak[0][0].theta <-10 && peak[0][0].snr > 3"); 

  cv->cd(5); 
  c.Draw("coherent[1][0].I", "peak[1][0].theta <0"); 
  cv->cd(6); 
  c.Draw("coherent[1][0].Q / coherent[1][0].I", "peak[1][0].theta <-10 && peak[1][0].snr > 3"); 
  cv->cd(7); 
  c.Draw("coherent[1][0].U / coherent[1][0].I", "peak[1][0].theta <-10 && peak[1][0].snr > 3"); 
  cv->cd(8); 
  c.Draw("coherent[1][0].V / coherent[1][0].I", "peak[1][0].theta <-10 && peak[1][0].snr > 3"); 




}
