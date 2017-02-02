{
gROOT->ProcessLine(".x macro/setupProof.C "); 
gROOT->ProcessLine(" .L macro/plotCuts.C "); 

TChain decimated("decimated"); decimated.Add("decimated/*.root"); decimated.SetProof();
TChain wais("wais"); wais.Add("wais/*.root"); wais.SetProof();

plotThermalCuts(&wais,&decimated,"cuts_deconvolved.pdf");

}

