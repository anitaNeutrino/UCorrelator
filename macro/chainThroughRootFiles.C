void chainThroughRootFiles(int data_start = 2, int data_stop=367) 
{
  //Step 1: load data
  TString tmp; 
  TChain* chain = new TChain("anita4"); 
  for (int i = data_start; i<= data_stop; i++){
    std::cout<< "load run "<< i <<std::endl;
    // tmp.Form("sparsedAllRuns/a4all/%d_max_999_sinsub_10_3_ad_2.root",i); 
    tmp.Form("/Volumes/SDCard/data/a4all/%d_max_30001_sinsub_10_3_ad_2.root",i); 
    chain->Add(tmp.Data());
  }
  // TFile *file = TFile::Open("sparsedAllRuns/sparsedAllRuns.root","RECREATE");
  TFile *file = TFile::Open("sparsedAllRuns/2-367_max_30001_sinsub_10_3_ad_2.root","RECREATE");
  chain->CloneTree(-1,"fast");
  file->Write();
}