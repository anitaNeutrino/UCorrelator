{
gROOT->ProcessLine(".x macro/setupProof.C "); 
gROOT->ProcessLine(" .L macro/plotCuts.C "); 

TChain a4all("anita4"); a4all.Add("/Volumes/SDCard/data/a4all/*30001*.root"); a4all.SetProof();
TChain wais("wais"); wais.Add("/Volumes/SDCard/data/wais/*30001*.root"); wais.SetProof();

plotThermalCuts(&wais,&a4all,"cuts_deconvolved.pdf");

}

