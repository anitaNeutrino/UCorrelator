
TH2 * eff2D(TTree * t, const char * var_x, int nbinsx, double xmin, double xmax, 
                       const char * var_y, int nbinsy, double ymin, double ymax, double scale, 
                       const char * cut = "weight * (O>=0)", const char * ordering = ">>") 
{
  TH2 * h  = new TH2D("2Ddist","MC Distribution", nbinsx, xmin, xmax, nbinsy, ymin, ymax); 
  TH2 * heff  = new TH2D("2Deff","MC Efficiency", nbinsx, xmin, xmax, nbinsy, ymin, ymax); 

  
  TString str;
  str.Form("%s:%s>>2Ddist",var_y,var_x); 
  int n = t->Draw(str.Data(), cut,"goff"); 



//  h->DrawCopy("colz"); 
 for (int i = 1; i <=h->GetNbinsX(); i++)
 {
   for (int j = 1; j <=h->GetNbinsY(); j++)
   {
        heff->SetBinContent(i,j, h->Integral(
                 ordering[0] == '>' ? i : 0, 
                 ordering[0] == '>' ? h->GetNbinsX()+1 : i, 
                 ordering[1] == '>' ? j : 0, 
                 ordering[1] == '>' ? h->GetNbinsY()+1 : j) / scale ); 
   }
 }

 heff->GetXaxis()->SetTitle(var_x); 
 heff->GetYaxis()->SetTitle(var_y); 
 delete h; 
 return heff; 
}
