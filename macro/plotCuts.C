#ifndef CUTS_LOADED
#include "macro/cuts.C"
#endif 



void plotThermalCuts(TChain * wais, 
                     TChain * decimated,
                     const char * output_file = 0, 
                     TH2 * tmpl = new TH2I("plot","Map Peak vs. Deconvolved Peak Hilbert; Deconvolved Peak Hilbert; Map Peak", 200, 0,0.05,100,0,.5), 
                     const char * draw= "peak.value:deconvolved.peakHilbert",
                     TCut wais_cuts = isReal && isWais, 
                     std::vector<TCut> thermalCuts = { isReal && !isWais ,brightestPeak, !isSun , !isNorth , aboveHorizon ,blastCut, badReconstructionCut , notTooFiltered,  triggered}, 
                     std::vector<const char*>  cut_labels = { "goodRF","brightest","!sun","!North", "upward",  " !blast", "goodReco", "!toofiltered", "triggered" } 
                     ) 

{


  int ncuts = thermalCuts.size(); 


  TCut cut = "" ;
  TString title = ""; 
  for (int i=0; i < ncuts; i++) 
  {
    cut += thermalCuts[i]; 
    TCanvas * c = new TCanvas(TString::Format("c_cut%d",i), cut,800,800); 
    TH2 * h = (TH2*) tmpl->Clone(TString::Format("%s_%d",tmpl->GetName(),i)); 
    TH2 * h2 = (TH2*) tmpl->Clone(TString::Format("wais_%s_%d",tmpl->GetName(),i)); 

    if (i > 0) title += ","; 
    if (cut_labels.size() > i)
    {
      title += cut_labels[i]; 
    }
    else
    {
      title += thermalCuts[i]; 
    }

    h->SetTitle(title); 
    c->cd(); 
    decimated->Draw(TString::Format("%s >> %s",draw, h->GetName()), cut, "colz"); 
    wais->Draw(TString::Format("%s >> %s",draw, h2->GetName()), wais_cuts, ""); 

    h->DrawCopy("colz"); 
    h2->SetMarkerColor(4); 
    h2->DrawCopy("same"); 

    if (output_file)
    {
      c->Print( TString::Format("%s%s", output_file, i == 0 ? "(" : i == ncuts -1 ? ")" : ""), "pdf");
    }

    delete h; 
    delete h2; 
  }


}

