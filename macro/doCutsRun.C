{

//  gROOT->ProcessLine(".x macro/setupProof.C "); 



  TChain wais("wais"); wais.Add("wais/*.root");

  for (int i = 130; i <=439; i++)
  {
    TChain decimated("decimated"); 
    decimated.Add(TString::Format("decimated/%d*.root",i));
    TString f = TString::Format("cutsByRun/cuts%d.pdf",i); 
    plotThermalCuts(&wais,&decimated,f.Data());
  }


}

