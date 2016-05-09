AnitaEventSummary * tryEvent(int run = 352, int event = 60832108, bool sinsub = false) 
{

//  AnalysisWaveform::enableDebug(true); 

  FFTtools::loadWisdom("wisdom.dat"); 

  gStyle->SetOptStat(0); 

  AnitaDataset d(run); 
  d.getEvent(event); 



  TFile out("test.root","RECREATE"); 
  FilterStrategy strategy(&out); 
  double fmins[]={0.23, 0.43}; 
  double fmaxs[]={0.29, 0.49}; 
  
  if (sinsub) 
  {
    strategy.addOperation(new UCorrelator::SineSubtractFilter(0.05, 0, 4)); 
  }
  else
  {
    UCorrelator::applyAbbysFilterStrategy(&strategy); 
  }

  printf("strategy applied!\n"); 

//  AnalysisWaveform::enableDebug(true); 
  printf("creating event\n"); 
  FilteredAnitaEvent * fae = new FilteredAnitaEvent(d.useful(), &strategy, d.gps(), d.header(),true); 



  printf("processed strategy!\n"); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  UCorrelator::Analyzer anal;
  anal.analyze(fae, sum); 





  /*
  TCanvas * all[2]; 
  all[0] = new TCanvas("horiz","horiz!",1800,1200); 
  all[1] = new TCanvas("vertic","vertic!",1800,1200); 

  all[0]->Divide(6,8); 
  all[1]->Divide(6,8); 
  */
  /*

  for (int i =0; i < 48; i++)
  {
    for (int pol = 0; pol < 2; pol++)
    {
      AnitaPol::AnitaPol_t p = (AnitaPol::AnitaPol_t)pol; 

      AnalysisWaveform * wf = fae->getFilteredGraph(i,p); 

  }

  printf("processed strategy!\n"); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  UCorrelator::Analyzer anal;
  anal.analyze(fae, sum); 





  /*
  TCanvas * all[2]; 
  all[0] = new TCanvas("horiz","horiz!",1800,1200); 
  all[1] = new TCanvas("vertic","vertic!",1800,1200); 

  all[0]->Divide(6,8); 
  all[1]->Divide(6,8); 
  */
  /*

  for (int i =0; i < 48; i++)
  {
    for (int pol = 0; pol < 2; pol++)
    {
      AnitaPol::AnitaPol_t p = (AnitaPol::AnitaPol_t)pol; 

      AnalysisWaveform * wf = fae->getFilteredGraph(i,p); 
      wf->even()->Write(TString::Format("g%d_%s", i, pol == 0 ? "hpol" : "vpol")); 
    }
  }
  */

//  all[0]->SaveAs("horiz.png"); 
//  all[1]->SaveAs("vertic.png"); 
  TCanvas * c1 = new TCanvas(); 
  c1->Divide(2,2); 


  UCorrelator::Correlator corr(180,0,360, 90, -60,25); 
  corr.setPadFactor(3); 

//  for (int i = 0; i < 10; i++) 
//  {
//    corr.compute(fae,AnitaPol::kVertical); 

   corr.compute(fae,AnitaPol::kHorizontal); 
//   corr.setGroupDelayFlag(0); 
   c1->cd(1); 
   corr.getHist()->DrawCopy("colz"); 

   UCorrelator::peakfinder::RoughMaximum rough[3]; 
   UCorrelator::peakfinder::FineMaximum fine[3]; 

   TCanvas * c2 = new TCanvas(); 
   c2->Divide(3,2); 
   int npeaks = UCorrelator::peakfinder::findIsolatedMaxima((const TH2D*)corr.getHist(), 20., 3, rough,true); 

   TH2D * zoomed = 0; 
   UCorrelator::WaveformCombiner comb(10); 
   for (int i = 0; i < npeaks; i++) 
   {
     c1->cd(2+i); 
     zoomed = corr.computeZoomed(rough[i].x, rough[i].y, 150, 0.2, 150, 0.2, 0,zoomed); zoomed->SetTitle(TString::Format("Zoomed peak %d", i)); zoomed->DrawCopy("colz"); 
     UCorrelator::peakfinder::doPeakFindingQuadratic25(zoomed, &fine[i]); 
     TMarker * m = new TMarker(fine[i].x, fine[i].y, 2); 
     double angle = 90. / TMath::Pi()* atan2(2*fine[i].covar, fine[i].sigma_x * fine[i].sigma_x - fine[i].sigma_y * fine[i].sigma_y);
     TEllipse *el = new TEllipse(fine[i].x, fine[i].y, fine[i].sigma_x, fine[i].sigma_y, 0, 360, angle); 
     el->SetFillStyle(0); 
     el->SetLineColor(3); 
     el->Draw(); 
     m->SetMarkerSize(2); 
     m->SetMarkerColor(3); 
     m->Draw(); 
     printf("--------------------------------------\n"); 
     printf("%f %f\n", fine[i].y, fine[i].x); 
     comb.combine(fine[i].x, -fine[i].y, fae, AnitaPol::kHorizontal); 
     c2->cd(1+i); 
     TString title; 
     title.Form("Coherent peak %d",i); 
     AnalysisWaveform * copy = new AnalysisWaveform(*comb.getCoherent()); copy->updateEven()->SetTitle(title);  copy->drawEven(); 
     c2->cd(4+i); copy->drawPowerdB();  
   }


// }
//  corr.getNorm()->DrawCopy("colz"); 
//  c1->SaveAs("norm.png"); 


//  corr.getHist()->DrawCopy("colz"); 
//  c1->SaveAs("corr.png"); 

//  */

//  sum.Print(); 

   /*

  TFile fcorr("fcorr.root","RECREATE"); 
  corr.getHist()->Write("corr"); 
  for (int i = 0 ; i < 48; i++) 
  {
    AnalysisWaveform * wf = fae->getFilteredGraph(i,AnitaPol::kHorizontal); 
    wf->even()->Write(TString::Format("g%d_%s", i, "hpol")); 
    wf = fae->getRawGraph(i,AnitaPol::kHorizontal); 
    wf->even()->Write(TString::Format("g%d_%s_unfiltered", i, "hpol")); 
    wf = fae->getFilteredGraphAtStage(i,AnitaPol::kHorizontal,0); 
    wf->even()->Write(TString::Format("g%d_%s_filtered_0", i, "hpol")); 
    wf = fae->getFilteredGraphAtStage(i,AnitaPol::kHorizontal,1); 
    wf->even()->Write(TString::Format("g%d_%s_filtered_1", i, "hpol")); 

    for (int j = i+1; j < 48; j++)
    {
      corr.getCorrelationGraph(i,j)->even()->Write(TString::Format("corr_%d_%d",i,j)); 
    }
  }
  */



  //corr.dumpDeltaTs("delta_ts.root");

  FFTtools::saveWisdom("wisdom.dat"); 
  return sum; 

}
