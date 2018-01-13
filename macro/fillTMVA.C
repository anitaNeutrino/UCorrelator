#include "AnitaTMVA.h" 

void _fillTMVA(const char * file, const char * tree="simulation", const char * weights="weights/thermal_cuts_Fisher.weights.xml", const char * name="F", const char * vars_file="tree_vars.tmva") 
{
  AnitaTMVA::MVAVarSet vars(vars_file); 

  TFile f(file,"update"); 
  TTree * t = (TTree*) f.Get(tree); 
  AnitaTMVA::evaluateTMVA(t, vars, name, weights); 
  t->Write(); 

}


void fillTMVA()
{

  for (int i = 1; i <=200; i+=40)
  {
    TString str; 
    str.Form("thermalTrees/simulated_%d-%d_max_1001_sinsub_10_3_ad_2.root", i, i + 39); 
    _fillTMVA(str.Data(),"simulation"); 
  }
return;
   for (int i = 120; i <=155; i+=40) 
  {
    TString str; 
    str.Form("thermalTrees/wais_%d-%d_max_30001_sinsub_10_3_ad_2.root", i, i + 39); 
    _fillTMVA(str.Data(),"wais"); 
  }

  for (int i = 50; i <=367; i+=40) 
  {
    TString str; 
    str.Form("thermalTrees/a4all_%d-%d_max_30001_sinsub_10_3_ad_2.root", i, i + 39); 
    _fillTMVA(str.Data(),"anita4"); 
  }


}
 

