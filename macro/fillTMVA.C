#include "AnitaTMVA.h" 

void fillTMVA(const char * file, const char * tree="simulation", const char * weights="thermal/weights/thermal_cuts_Fisher.weights.xml", const char * name="F", const char * vars_file="tree_vars.tmva") 
{
  AnitaTMVA::MVAVarSet vars(vars_file); 

  TFile f(file,"update"); 
  TTree * t = (TTree*) f.Get(tree); 
  AnitaTMVA::evaluateTMVA(t, vars, name, weights); 
  t->Write(); 

}


void fillAllTMVA()
{

  for (int i = 1; i <=401; i+=100) 
  {
    TString str; 
    str.Form("thermalTrees/simulated_%d-%d_sinsub_10_3_ad_2.root", i, i + 99); 
    fillTMVA(str.Data(),"simulation"); 
  }

  for (int i = 130; i <=430; i+=10) 
  {
    TString str; 
    str.Form("thermalTrees/a3all_%d-%d_sinsub_10_3_ad_2.root", i, i + 9); 
    fillTMVA(str.Data(),"anita3"); 
  }


}
 

