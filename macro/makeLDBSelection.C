
void makeLDBSelection(int start_run=130, int end_run=165, int nsbins = 1e4, int nsecs = 60, int thresh=20 )
{ 
  TChain c("headTree"); 

  for (int i = start_run; i<= end_run; i++) 
  {
    c.Add(TString::Format("%s/run%d/headFile%d.root", getenv("ANITA_ROOT_DATA"),i,i)); 
  }

  double t0 = c.GetMinimum("triggerTime"); 
  double t1 = c.GetMaximum("triggerTime")+1; 

  int ntimebins = ceil((t1-t0)/nsecs); 

  TFile f("data/ldbSelection.root","RECREATE"); 

  TH2D *ldbHist = new TH2D("ldbHist","LDBHist",  ntimebins, t0, t0 + ntimebins * nsecs, nsbins, 0,1e9); 


  c.Draw("triggerTimeNs:triggerTime >> ldbHist","trigType &1","goff"); 

  ldbHist->DrawCopy("colz"); 
  printf("Killing all bins below %d\n",thresh); 


  for (int i = 1; i <= ldbHist->GetNbinsX(); i++)
  {
    for (int j = 1; j <= ldbHist->GetNbinsY(); j++)
    {
      if (ldbHist->GetBinContent(i,j) < thresh)
      {
        ldbHist->SetBinContent(i,j,0); 
      }
    }
  }


  f.Write(); 

}
