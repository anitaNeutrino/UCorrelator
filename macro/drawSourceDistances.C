const int nkeys =6; 
const char * keys[nkeys] = {"txs","fava","sn","grb","grb_12","grb_24"}; 

void drawSourceDistances(const char * only=0)
{

  TFile f("source_distances.root"); 

  for (int ikey = 0; ikey < nkeys; ikey++)
  {
    if (only && !strstr(only,keys[ikey]))
      continue; 

    TH3 * h = (TH3*) f.Get(keys[ikey]); 
    TH3 * h_raw = (TH3*) f.Get(Form("%s_raw",keys[ikey])); 
    TH3 * h_int = (TH3*) f.Get(Form("%s_int",keys[ikey])); 
    printf ("%p %p %p\n", h, h_raw, h_int); 

    int nz = h->GetNbinsZ(); 


    TCanvas * c = new TCanvas(Form("c%s", keys[ikey]), keys[ikey], 1800,900); 
    TCanvas * c_int = new TCanvas(Form("c%s_int", keys[ikey]), Form("Interpolated %s",keys[ikey]), 1800,900); 
    TCanvas * c_raw = new TCanvas(Form("c%s_raw", keys[ikey]), Form("Raw %s", keys[ikey]), 1800,900); 

    printf ("%p %p %p\n", c, c_raw, c_int); 

    if (nz > 4) 
    {
      c->Divide(4,ceil(nz/4.)); 
      c_int->Divide(4,ceil(nz/4.)); 
      c_raw->Divide(4,ceil(nz/4.)); 
    }
    else if (nz > 1) 
    {
      c->Divide(2,ceil(nz/2.)); 
      c_int->Divide(2,ceil(nz/2.)); 
      c_raw->Divide(2,ceil(nz/2.)); 
    }

    for (int iz = 1; iz <= nz; iz++) 
    {
      printf ("%d\n",iz); 

      h->GetZaxis()->SetRange(iz,iz); 

      h_int->GetZaxis()->SetRange(iz,iz); 

      h_raw->GetZaxis()->SetRange(iz,iz); 

      c->cd(nz > 1 ? iz : 0); 
      auto proj =  h->Project3D("yx"); 
      proj->SetName(Form("proj_%s",keys[ikey]));
      auto proj_raw =  h_raw->Project3D("yx"); 
      proj_raw->SetName(Form("proj_raw_%s",keys[ikey]));
      auto proj_int =  h_int->Project3D("yx"); 
      proj_int->SetName(Form("proj_int_%s",keys[ikey]));
      printf("%p %p %p\n", proj, proj_raw, proj_int); 

      TString title = Form("Slice %d , t #in [%f %f] ;RA;dec", iz, h->GetZaxis()->GetBinLowEdge(iz), h->GetZaxis()->GetBinLowEdge(iz+1));
      printf("%s\n",title.Data()) ;
      proj->SetTitle(title); 
      proj->SetStats(0); 
      proj->DrawCopy("col2z"); 



      proj_int->SetTitle(title); 
      proj_int->SetStats(0); 
      c_int->cd(nz > 1 ? iz : 0); 
      proj_int->DrawCopy("col2z"); 



      proj_raw->SetTitle(title); 
      proj_raw->SetStats(0); 

      c_raw->cd(nz > 1 ? iz : 0); 
      proj_raw->DrawCopy("colz"); 
    }
  }



}
