// root plot trend of singlets.C
// The output will be in plots_trendOfSinglets.root
void plotTrendOfSinglets(const char * var = "deconvImpulsivity")
{
  TChain chain("trend"); 
  chain.Add("source_maps/anita4/trendOfSinglets_3.5sigma_30002_mod30_50_367.root");
  TFile plots_trendOfSinglets("plots_trendOfSinglets.root","RECREATE");


  chain.Draw("percentOfData:N_singlets/percentOfData:N_singlets_nearbase/percentOfData:N_singlets_notnearbase/percentOfData");
  TGraph *gr1 = new TGraph(chain.GetSelectedRows(), chain.GetV1(), chain.GetV2());
  TGraph *gr2 = new TGraph(chain.GetSelectedRows(), chain.GetV1(), chain.GetV3());
  TGraph *gr3 = new TGraph(chain.GetSelectedRows(), chain.GetV1(), chain.GetV4());
  gr1->SetName("gr1");
  gr2->SetName("gr2");
  gr3->SetName("gr3");
  gr1->GetXaxis()->SetRangeUser(0.,1.1);
  gr1->GetYaxis()->SetRangeUser(0.,150);
  gr1->SetMarkerStyle(3);
  gr2->SetMarkerStyle(3);
  gr3->SetMarkerStyle(3);
  gr1->SetMarkerColor(1);
  gr2->SetMarkerColor(2);
  gr3->SetMarkerColor(3);
  TCanvas * dists = new TCanvas("c1","c1"); 

  gr1->Draw("ap"); //draw graph in current pad
  gr2->Draw("p"); //draw graph in current pad
  gr3->Draw("p"); //draw graph in current pad
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry("gr1","N_singlets","p");
  legend->AddEntry("gr2","N_singlets_nearbase","p");
  legend->AddEntry("gr3","N_singlets_notnearbase","p");
  legend->Draw();
  gr1->GetXaxis()->SetTitle("ratio of data");
  gr1->GetYaxis()->SetTitle("N of singlets");
  plots_trendOfSinglets.cd(); 
  gr1->Write(); 
  gr2->Write(); 
  gr3->Write(); 
  legend->Write();
  std::cout << "finished"<<std::endl;

}
