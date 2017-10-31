#include "TFile.h"
#include "TTree.h"
#include "AnitaEventSummary.h" 
#include "FFTtools.h"
// plot delta theta vs phi for all wais event runs.
void plotMultiSummary()
{
	TString PlotPrefix = "Plot_";
	TString pointDir = "../drivers/wais/";
	TString fEnd = "_max_30007_sinsub_10_3_ad_2.root";	
	// TString fEnd = "_max_30002_.root";	
	TString outf = pointDir + PlotPrefix + fEnd;
	
	TChain* chain = new TChain("wais");
	for(int i = 123; i < 153; i++)
	// for(int i = 139; i < 141; i++)
	{
		// if(i==129 || i==130 || i==131 || i==132) continue;
		//if(i==132) continue;
		TString inf = pointDir + TString::Itoa(i, 10) + fEnd;
		chain->Add(inf);
	}

	TFile *file = TFile::Open("../drivers/wais/allWais" + fEnd,"RECREATE");
	  chain->CloneTree(-1,"fast");
	  file->Write();
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
	TH3D* dTheta_dPhi_Power_H = new TH3D("dTheta_dPhi_Power_H", "dTheta_dPhi_Power_H", 100, -180, 180, 100, -25, 5, 100, 0, 0.15);	
	TH3D* dTheta_dPhi_Power_V = new TH3D("dTheta_dPhi_Power_V", "dTheta_dPhi_Power_V", 100, -180, 180, 100, -25, 5, 100, 0, 0.15);	
	// TH3D* dTheta_dPhi_RMS_H = new TH3D("dTheta_dPhi_RMS_H", "dTheta_dPhi_RMS_H", 100, -180, 180, 100, -25, 5, 100, 0, 50);	
	// TH3D* dTheta_dPhi_RMS_V = new TH3D("dTheta_dPhi_RMS_V", "dTheta_dPhi_RMS_V", 100, -180, 180, 100, -25, 5, 100, 0, 50);	
	TH2D* dThetavsdPhiH = new TH2D("dThetavsdPhiH", "dThetavsdPhiH", 100, -2, 3, 100, -2, 2);	
	TH2D* dThetavsdPhiV = new TH2D("dThetavsdPhiV", "dThetavsdPhiV", 100, -2, 3, 100, -2, 2);	

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
		// double theta1 = sum->mostImpulsivePeak().theta;
		// double phi1 = sum->mostImpulsivePeak().phi;
		double dist = w.distance/1e3;
		// double dtheta = theta1 - waistheta;
		// double dphi = FFTtools::wrap(phi1 - waisphi, 360, 180);
		
		if(pulser == 1){
			// if(fabs(dtheta) > 3 or fabs(dphi)>5 ){
			// 	continue;
			// }
			// dThetavPhiH->Fill(waisphi, dtheta);
			dThetavPhiH->Fill(sum->wais.phi, sum->peak[0][0].theta - sum->wais.theta);
			dThetavPhiAll->Fill(sum->wais.phi, sum->peak[0][0].theta - sum->wais.theta);
			for(int ant=0; ant < 48 ; ant++){
				dTheta_dPhi_Power_H->Fill(FFTtools::wrap(sum->wais.phi - (sum->channels[0][ant].getPhi() - 2)*22.5,360,0), 10 - sum->wais.theta , sum->channels[0][ant].avgPower);
				// dTheta_dPhi_RMS_H->Fill(FFTtools::wrap(sum->wais.phi - (sum->channels[0][ant].getPhi() - 2)*22.5,360,0), 10 - sum->wais.theta , sum->channels[0][ant].rms);
			}

			dThetavsdPhiH->Fill(FFTtools::wrap(sum->peak[0][0].phi - sum->wais.phi, 360, 0), sum->peak[0][0].theta - sum->wais.theta);
		}
		if(pulser == 5){
			dThetavPhiV->Fill(sum->wais.phi, sum->peak[1][0].theta - sum->wais.theta);
			dThetavPhiAll->Fill(sum->wais.phi, sum->peak[1][0].theta - sum->wais.theta);
			// dThetavPhiV->Fill(waisphi, dtheta);
			for(int ant=0; ant < 48 ; ant++){
				dTheta_dPhi_Power_V->Fill(FFTtools::wrap(sum->wais.phi - (sum->channels[1][ant].getPhi() - 2)*22.5,360,0), 10 - sum->wais.theta , sum->channels[1][ant].avgPower);
				// dTheta_dPhi_RMS_V->Fill(FFTtools::wrap(sum->wais.phi - (sum->channels[1][ant].getPhi() - 2)*22.5,360,0), 10 - sum->wais.theta , sum->channels[1][ant].rms);
			}
			dThetavsdPhiV->Fill(FFTtools::wrap(sum->peak[1][0].phi - sum->wais.phi, 360, 0), sum->peak[1][0].theta - sum->wais.theta);
		}
		dThetavPhiAll->GetXaxis()->SetTitle("wais.phi/[degree]");
		dThetavPhiAll->GetYaxis()->SetTitle("deltaTheta/[degree]");

		dThetavsdPhiH->GetXaxis()->SetTitle("dPhi/[degree]");
		dThetavsdPhiH->GetYaxis()->SetTitle("dTheta/[degree]");

		dThetavsdPhiV->GetXaxis()->SetTitle("dPhi/[degree]");
		dThetavsdPhiV->GetYaxis()->SetTitle("dTheta/[degree]");

		
	}
	TFile f2(outf.Data(), "RECREATE");
	f2.cd();
	dThetavPhiH->Write();
	dThetavPhiV->Write(); 
	dThetavPhiAll->Write();
	dTheta_dPhi_Power_H->Write();
	dTheta_dPhi_Power_V->Write();
	// dTheta_dPhi_RMS_H->Write();
	// dTheta_dPhi_RMS_V->Write();
	dThetavsdPhiH->Write();
	dThetavsdPhiV->Write();

	dTheta_dPhi_Power_H->FitSlicesZ();
	dTheta_dPhi_Power_V->FitSlicesZ();
	// dTheta_dPhi_RMS_H->FitSlicesZ();
	// dTheta_dPhi_RMS_V->FitSlicesZ();

	gDirectory->Get("dTheta_dPhi_Power_H_1")->Write();
	gDirectory->Get("dTheta_dPhi_Power_V_1")->Write();
	// gDirectory->Get("dTheta_dPhi_RMS_H_1")->Write();
	// gDirectory->Get("dTheta_dPhi_RMS_V_1")->Write();


	f2.Close();
}
