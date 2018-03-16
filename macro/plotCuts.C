#ifndef CUTS_LOADED
#include "cuts.C"
#endif 

void plotThermalCuts(TChain * chain, const char * output_file = 0, std::vector<TCut> cutsList ={}, std::vector<const char*>  cut_labels={}, bool isWais=0,bool isMC =0, bool isSequential = 0){
  std::cout<<"############################################"<<std::endl;
  int ncuts = cutsList.size(); 
  TCut weight = "mc.weight";
  TCut cut,cut0;
  TString title;
  if(isWais){
    //initial cut for wais
    cut=allWais;
    cut0=allWais;
    title="allWais";
  }else{
     // intitial cut for other events.
    cut="1";
    cut0="1";
    title="1";
  }
  double count0 = 0;
  double frac = 0;
  for (int i=-1; i < ncuts; i++) 
  {
    if (i != -1){
      if(isSequential){
        cut += cutsList[i];
        if (i>=0) title += ","; 
        if (cut_labels.size() > i)
        {
          title += cut_labels[i]; 
        }
        else
        {
          title += cutsList[i]; 
        } 
      }else{
        cut = cut0 + cutsList[i];
        title = cut_labels[i];
      }
    }

    TCanvas * canvas = new TCanvas(TString::Format("c_cut%d",i), cut,800,800); 
    // TH2D is double typed plot then we can count the weighted event.
    TH2 * h = new TH2D("plot","Vpol impulsivity:Hpol impulsivity; Hpol impulsivity; Vpol impulsivity;", 300, 0,1,300,0,1); 
    h->SetTitle(title); 
    canvas->cd();
    double count;
    if(isMC){
      count = chain->Draw(TString::Format("%s >> %s","deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure", h->GetName()),cut*weight, "colz");
      count = h->Integral();
    }else{
      count = chain->Draw(TString::Format("%s >> %s","deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure", h->GetName()),cut, "colz");
    } 
    if(i == -1){
      count0 = count;
    } 
    std::cout<< title<<" \t"<< count<<"("<< 100*double(count)/double(count0)<<"%)"<<std::endl;
    // h->SetMarkerColor(4); 
    // h->DrawCopy("colz"); 

    if (output_file)
    {
      canvas->Print( TString::Format("%s%s", output_file, i == -1 ? "(" : i == ncuts -1 ? ")" : ""), "pdf");
    }

    delete h; 
  }
}

void plotCuts(){
  // bool isSequential or if only has one cut each time
  bool isSequential=0;
  std::vector<TCut> cutsList;
  std::vector<const char*>  cut_labels;

  // TChain wais("wais"); wais.Add("/Volumes/SDCard/data/wais/*30001*.root");
  // cutsList = {isReal, notGlitch, notBadReconstruction, notBlast, triggered, notMasked,notStrongCW,notHical,goodPointingWais};
  // cut_labels = {"isReal", "notGlitch", "notBadReconstruction", "notBlast", "triggered", "notMasked","notStrongCW","notHical","goodPointingWais"};
  // plotThermalCuts(&wais,"cuts_deconvolved_wais.pdf",cutsList,cut_labels,1,0,isSequential);


  TChain a4all("anita4"); a4all.Add("/Volumes/SDCard/data/a4all/*30002*.root");
  cutsList = {isReal, notGlitch, notBadReconstruction, notBlast, triggered, notMasked,notStrongCW,notHical};
  cut_labels = {"isReal", "notGlitch", "notBadReconstruction", "notBlast", "triggered", "notMasked","notStrongCW","notHical"};
  plotThermalCuts(&a4all,"cuts_deconvolved_a4all.pdf",cutsList,cut_labels,0,0,isSequential);

  // // TChain mc("simulation"); mc.Add("/Volumes/SDCard/data/simulated/*1000*.root"); // Energy 222
  // TChain mc("simulation"); mc.Add("/Volumes/SDCard/data/simulated/*1001*.root"); // Min Bias Eneryg 222
  // cutsList = {isReal, notGlitch, notBadReconstruction, notBlast, triggered, notMasked,notStrongCW, notHical,goodPointingMC};
  // cut_labels = {"isReal", "notGlitch", "notBadReconstruction", "notBlast", "triggered", "notMasked","notStrongCW", "notHical","goodPointingMC"};
  
  // // // cutsList = {isReal,goodPointingMC};
  // // // cut_labels = {"isReal", "goodPointingMC"};
  // plotThermalCuts(&mc,"cuts_deconvolved_simulated.pdf",cutsList,cut_labels,0,1,isSequential);
}



