#include "TFile.h"
#include "TTree.h"
#include "AnitaEventSummary.h" 
#include "FFTtools.h"
// plot delta theta vs phi for all wais event runs.
void plotMultiSummary()
{
	TString PlotPrefix = "Plot_";
	TString pointDir = "../drivers/a4all/";
	TString fEnd = "_max_30000_sinsub_10_3_ad_2.root";	
	TString outf = pointDir + PlotPrefix + fEnd;
	
	TChain* chain = new TChain("anita4");
	for(int i = 123; i < 153; i++)
	// for(int i = 139; i < 141; i++)
	{
		if(i==129 || i==130 || i==131 || i==132) continue;
		//if(i==132) continue;
		TString inf = pointDir + TString::Itoa(i, 10) + fEnd;
		chain->Add(inf);
	}

	AnitaEventSummary* sum = new AnitaEventSummary;
	AnitaEventSummary::PointingHypothesis point;
	AnitaEventSummary::SourceHypothesis w;
	chain->SetBranchAddress("summary", &sum);

	Long_t nEntries = chain->GetEntries();
	chain->GetEntry(0);
	Double_t firstEvent = sum->eventNumber;
	chain->GetEntry(nEntries-1);
	Double_t lastEvent = sum->eventNumber;
	TH2D* dThetavPhiH = new TH2D("HthetaPhi", "HthetaPhi", 200, 0, 360, 100, -2, 2);	
	TH2D* dThetavPhiV = new TH2D("VthetaPhi", "VthetaPhi", 200, 0, 360, 100, -2, 2);	
	TH2D* dThetavPhiAll = new TH2D("thetaPhi", "thetaPhi", 200, 0, 360, 100, -2, 2);	

	for (int j = 0; j < nEntries; j ++)
	{
		chain->GetEntry(j);
		int pulser = sum->flags.pulser;
		if(pulser == 1) {
			point = sum->peak[0][0];
		}else if(pulser == 5){
			point = sum->peak[1][0];
		}else{
			continue;
		}
		w = sum->wais;
		double waistheta = w.theta;
		double waisphi = w.phi;
		double theta1 = point.theta;
		double phi1 = point.phi;
		double dist = w.distance/1e3;
		double dtheta = theta1 - waistheta;
		double dphi = FFTtools::wrap(phi1 - waisphi, 360, 180);
		if(fabs(dtheta) > 2 ){
			continue;
		}
		if(pulser == 1){
			dThetavPhiH->Fill(waisphi, dtheta);
		}
		if(pulser == 5){
			dThetavPhiV->Fill(waisphi, dtheta);
		}
		dThetavPhiAll->Fill(waisphi, dtheta);
		dThetavPhiAll->GetXaxis()->SetTitle("wais.phi/[degree]");
		dThetavPhiAll->GetYaxis()->SetTitle("deltaTheta/[degree]");

		
	}
	TFile f2(outf.Data(), "RECREATE");
	f2.cd();
	dThetavPhiH->Write();
	dThetavPhiV->Write(); 
	dThetavPhiAll->Write();
	f2.Close();
}
