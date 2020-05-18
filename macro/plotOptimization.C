

void plotOptimization(const char * key = "txs", const char* thermalAddString = "thermalTrees/simulated_txs_src_1-100*.root", int which_prior=1, int which_cl = 0)
{
  TFile f(Form("%s_opt.root", key)); 

  TH3 * eff = (TH3*) f.Get(Form("%s__eff", key)); 
  TH3 * bg = (TH3*) f.Get(Form("%s__bg", key)); 

  const int ncls =2;
  const int npriors = 2;
  double CLs[ncls] = {0.9,0.975}; 
  double priors[ncls] = {0.5,1}; 

  TH3 * sups[ncls][npriors]; 
  double Ds[ncls][npriors];
  double Fs[ncls][npriors];
  double logOs[ncls][npriors];

  for (int iprior = 0; iprior < npriors; iprior++) 
  {
    if (which_prior >=0 && iprior != which_prior) continue; 

    for (int icl = 0; icl < ncls; icl++) 
    {
      if (which_cl >=0 && icl != which_cl) continue; 


      double CL = CLs[icl]; 
      double prior = priors[iprior]; 

      TH3 * sup = (TH3*) f.Get(Form("%s_scaled_sup_%g_%g", key, CL*1000, prior*10)); 

      int minbin = sup->GetMinimumBin(); 
      int bx,by,bz; 
      sup->GetBinXYZ(minbin,bx,by,bz); 
      
      double better_threshold = 0.005; 
      double sup_estimate = sup->GetBinContent(bx,by,bz); 

      while (true) 
      {
        double sup_prime = sup->GetBinContent(bx-1,by,bz); 
        if ( (sup_prime-sup_estimate)/sup_estimate < better_threshold) 
        {
          bx--; 
          sup_estimate = sup->GetBinContent(bx,by,bz); 
          continue;
        }

        sup_prime = sup->GetBinContent(bx,by-1,bz); 
        if ((sup_prime-sup_estimate)/sup_estimate < better_threshold) 
        {
          by--; 
          sup_estimate = sup->GetBinContent(bx,by,bz); 
          continue;
        }

        if (bz == 100) break; 
        sup_prime = sup->GetBinContent(bx,by,bz+1); 
        if ((sup_prime-sup_estimate)/sup_estimate < better_threshold) 
        {
          bz++; 
          sup_estimate = sup->GetBinContent(bx,by,bz); 
          continue;
        }
        break; 
      }





      double eff_estimate = eff->GetBinContent(bx,by,bz); 
      double bg_estimate = bg->GetBinContent(bx,by,bz); 
      double F = eff->GetXaxis()->GetBinLowEdge(bx);
      double logO = eff->GetYaxis()->GetBinLowEdge(by);
      double D = eff->GetZaxis()->GetBinLowEdge(bz);

      Fs[icl][iprior] = F; 
      Ds[icl][iprior] = D; 
      logOs[icl][iprior] = logO; 

      int bins[] = {bz};
      int nbins =sizeof(bins)/sizeof(int); 

      printf("================\n"); 
      printf("Prior=%g\n", prior);
      printf("CL=%g\n", CL);
      printf("Best cut:\n"); 
      printf("\tF < %f\n", F);
      printf("\t-log10(O) < %f\n", logO);
      printf("\tD < %f\n", D);
      printf("Efficiency estimate: %f%%\n", 100*eff_estimate);
      int nobs = 100*bg_estimate; 
      int ntotal = 100; 

      double bg_up = ROOT::Math::gamma_quantile(0.84, nobs*7./100+prior, 1./7); 
      double bg_mid = ROOT::Math::gamma_quantile(0.5, nobs*7./100+prior, 1./7); 
      double bg_dn = ROOT::Math::gamma_quantile(0.16, nobs*7./100+prior, 1./7); 

      printf("Background estimate: %f (+%f, -%f) \n", bg_mid, (bg_up-bg_mid), bg_mid-bg_dn );
      printf("S_{up}(%f%%, %s prior): %f\n", CL*100, prior == 0.5 ? "Jeffreys" : "Uniform", sup_estimate); 

      for (int ibin = 0; ibin < nbins; ibin++) 
      {
        TMarker * m = 0;
        if (ibin == 0) 
        {
          m = new TMarker(F,logO,5); 
          m->SetMarkerSize(3); 
          m->SetMarkerColor(3); 

        }

        int bin = bins[ibin]; 
        D = eff->GetZaxis()->GetBinLowEdge(bin); 
        auto c = new TCanvas(Form("cbin_%d_%d_%d",ibin,icl,iprior),Form("D<=%f, prior=%g, CL=%g", D,prior,CL), 1200, 500);  ; 
        c->Divide(3,1); 
        c->cd(1); 
        eff->GetZaxis()->SetRange(bin,bin); 
        auto proj_eff = eff->Project3D("yx"); 
        proj_eff->SetTitle(Form("%s efficiency, D < %g;F;-logO", key, D)); 
        proj_eff->SetStats(false);
//        proj_eff->SetMaximum(1); 
//        proj_eff->SetMinimum(0); 
        proj_eff->DrawCopy("col2z"); 
        if (m) m->Draw(); 

        c->cd(2); 
        bg->GetZaxis()->SetRange(bin,bin); 
        auto proj_bg = bg->Project3D("yx"); 
        proj_bg->SetTitle(Form("%s bg, D < %g;F;-logO", key, D)); 
        proj_bg->SetStats(false);
        proj_bg->DrawCopy("col2z"); 
        if (m) m->Draw(); 
        gPad->SetLogz(); 

        c->cd(3); 
        sup->GetZaxis()->SetRange(bin,bin); 
        auto proj_sup = sup->Project3D("yx"); 
        proj_sup->SetTitle(Form("%s s_{up} (%f%%, %s prior), D < %g;F;-logO", key, CL, prior == 0.5 ? "Jeffreys" : "Uniform",  D)); 
        proj_sup->SetMaximum(proj_sup->GetMinimum()*1.2); 
        proj_sup->SetStats(false);
        proj_sup->DrawCopy("col2z"); 
        if (m) m->Draw(); 
      }
    }
  }


 

  if (thermalAddString) 
  {

    //now make Efficiency
    TCanvas * ceff = new TCanvas("eff","Efficiency"); 

    TChain toverlap("overlap"); 
    toverlap.Add(Form("%s_overlaps.root",key)); 

    TChain tthermal("simulation"); 
    tthermal.Add(thermalAddString); 

    TMultiGraph * mgE = new TMultiGraph; 
    TMultiGraph * mgSNR = new TMultiGraph; 
    mgE->SetTitle("Efficiency vs. Energy (binned, not fixed!)" ); 
    mgE->GetXaxis()->SetTitle("log10 (E /eV)"); 
    mgE->GetYaxis()->SetTitle("Efficiency"); 

    mgSNR->SetTitle("Efficiency vs. Coherent SNR (not including peak-finding inefficiency!)" ); 
    mgSNR->GetXaxis()->SetTitle("Coherent SNR"); 
    mgSNR->GetYaxis()->SetTitle("Efficiency "); 
 
    int igraphs = 0; 

    for (int icl = 0; icl < ncls; icl++) 
    {
       if (which_cl >=0 && icl != which_cl) continue; 
      for (int iprior = 0; iprior < npriors; iprior++) 
      {
        if (which_prior >=0 && iprior != which_prior) continue; 

        TH1D denom_E("denom_E","E Denominator",7,17.75,21.25);
        TH1D num_E("num_E","E Numerator",7,17.75,21.25);

        double F = Fs[icl][iprior]; 
        double logO = logOs[icl][iprior]; 
        double D = Ds[icl][iprior]; 
        
        //evaluate efficiency only with EVEN events 
        tthermal.Draw("log10(energy) >> denom_E","weight * isMostImpulsive * ((eventNumber_nnn & 1) == 0)","goff"); 
        toverlap.Draw("log10(mcE) >> num_E" ,Form("weight * ((event & 1) == 0) * (F > %g ) * (D_%s < %f) * ( O < %f)", F,key, D, pow(10,-logO)),"goff"); 

        TH1D denom_SNR("denom_snr","SNR Denominator",10,0,50); 
        TH1D num_SNR("num_snr","SNR Numerator",10,0,50); 

        tthermal.Draw("coherentSNR >> denom_snr","weight * isMC * ((eventNumber_nnn & 1) == 0)", "goff"); 
        toverlap.Draw("cohSNR>> num_snr" ,Form("weight * ((event & 1) == 0) * (F > %g ) * (D_%s < %f) * ( O < %f)", F,key, D, pow(10,-logO)),"goff"); 


        TGraphAsymmErrors *gE = new TGraphAsymmErrors; 
        TGraphAsymmErrors *gSNR = new TGraphAsymmErrors; 

        gE->BayesDivide(&num_E,&denom_E);
        gSNR->BayesDivide(&num_SNR,&denom_SNR);

        gE->SetTitle(Form("CL=%g, prior = %g", CLs[icl], priors[iprior])); 
        mgE->Add(gE); 

        gSNR->SetTitle(Form("CL=%g, prior = %g", CLs[icl], priors[iprior])); 
        mgSNR->Add(gSNR); 
        igraphs++; 
      }

    }

    ceff->Divide(2,1); 

    ceff->cd(1); 
    mgE->Draw("a plc pmc"); 
    mgE->GetYaxis()->SetLimits(0,1);
    if (igraphs > 1) 
      gPad->BuildLegend(0.4,0.1,0.9,0.4); 
    ceff->cd(2); 
    mgSNR->GetYaxis()->SetLimits(0,1);
    mgSNR->Draw("a plc pmc"); 
    if (igraphs > 1) 
      gPad->BuildLegend(0.4,0.1,0.9,0.4); 

  }


}
