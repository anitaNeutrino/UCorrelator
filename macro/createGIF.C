#include "TFile.h"
#include "TTree.h"
#include "AnitaEventSummary.h" 
// Given summary file, it can generate a GIF format movie for summary->peak[][] theta vs phi.  
//Author: Peng Cao
//root createGIF.root
void createGIF() 
{
   gSystem->Unlink("../drivers/all/movie.gif"); // delete old file
   TCanvas* c = new TCanvas();
   TFile *f = new TFile("../drivers/all/140_max_10000_deconv.root");
   if (f == 0) {
      printf("Error: cannot open root file.\n");
      return;
   }
   TTree *anita3 = (TTree*) f->Get("anita3");
   AnitaEventSummary *summary = NULL;     
   anita3->SetBranchAddress("summary",&summary); 
   int num_entries = anita3->GetEntries();
   cout << "number of entries is " << num_entries << endl; ;

   TH2D* hist= NULL;
   anita3->GetEntry(0);
   int startEventNumber= summary->eventNumber;
   char* histName = new char[20];
   int step = 10000;
   for (int ientry = 0; ientry < num_entries; ientry+=step) 
   {
         gStyle->SetOptStat(kFALSE);
         anita3->GetEntry(ientry); // after doing get entry you can forget about the tree 
         sprintf(histName,"%d",summary->realTime);
         hist = new TH2D("h2",histName,360,0.0,360.0,300,-90.0,90.0);
         anita3->Draw("-1*summary->peak[0][].theta:(summary->peak[0][].phi- pat.heading +360)%360>>h2","","colz",step,ientry);
         hist->GetXaxis()->SetTitle("phi");
         hist->GetYaxis()->SetTitle("theta");
         c->Print("../drivers/all/movie.gif+");
         delete hist;
       
      
   }
   c->Print("../drivers/all/movie.gif++");
   delete c;
   printf("movie.C ends\n"); 
   gApplication->Terminate(0); 
}