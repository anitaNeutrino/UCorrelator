#include "TFile.h"
#include "TTree.h"
#include "AnitaEventSummary.h" 
#include "FFTtools.h"
// plot delta theta vs phi for all wais event runs.
void _plotMultiSummary(TString treeName, TString pointDir, TString fEnd ,TString outf, int start_run, int end_run)
{	
	TChain* chain = new TChain(treeName);
	for(int i = start_run; i < end_run; i++)
	// for(int i = 139; i < 141; i++)
	{
		if (treeName.EqualTo("wais")){
			if(i==129 || i==130 || i==131 || i==132) continue;
		};
		TString inf = pointDir + TString::Itoa(i, 10) + fEnd;
		chain->Add(inf);
	}

	// TFile *file = TFile::Open("../drivers/wais/allWais" + fEnd,"RECREATE");
	//   chain->CloneTree(-1,"fast");
	//   file->Write();
	AnitaEventSummary* sum = new AnitaEventSummary;
	AnitaEventSummary::SourceHypothesis w;
	chain->SetBranchAddress("summary", &sum);

	Long_t nEntries = chain->GetEntries();
	// TH2D* dThetavPhiH = new TH2D("HthetaPhi", "HthetaPhi", 200, 0, 360, 100, -2, 2);	
	// TH2D* dThetavPhiV = new TH2D("VthetaPhi", "VthetaPhi", 200, 0, 360, 100, -2, 2);	
	// TH2D* dThetavPhiAll = new TH2D("thetaPhi", "thetaPhi", 200, 0, 360, 100, -2, 2);
	// TH3D* dTheta_dPhi_Power_H = new TH3D("dTheta_dPhi_Power_H", "dTheta_dPhi_Power_H", 100, -180, 180, 100, -25, 5, 100, 0, 0.15);	
	// TH3D* dTheta_dPhi_Power_V = new TH3D("dTheta_dPhi_Power_V", "dTheta_dPhi_Power_V", 100, -180, 180, 100, -25, 5, 100, 0, 0.15);	
	// TH3D* dTheta_dPhi_RMS_H = new TH3D("dTheta_dPhi_RMS_H", "dTheta_dPhi_RMS_H", 100, -180, 180, 100, -25, 5, 100, 0, 50);	
	// TH3D* dTheta_dPhi_RMS_V = new TH3D("dTheta_dPhi_RMS_V", "dTheta_dPhi_RMS_V", 100, -180, 180, 100, -25, 5, 100, 0, 50);	
	TH2D* dThetavsdPhiH = new TH2D("dThetavsdPhiH", "dThetavsdPhiH", 100, -5, 5, 250, -5, 5);	
	TH2D* dThetavsdPhiV = new TH2D("dThetavsdPhiV", "dThetavsdPhiV", 100, -5, 5, 250, -5, 5);	
	TH2D* dThetavsdPhi = new TH2D("dThetavsdPhi", "dThetavsdPhi", 100, -5, 5, 250, -5, 5);	
	TH2D* dThetavsSNR = new TH2D("dThetavsSNR", "dThetavsSNR", 30, 0, 100, 50, -2, 2);	
	TH2D* dPhivsSNR = new TH2D("dPhivsSNR", "dPhivsSNR", 30, 0, 100, 50, -5, 5);	

	for (int j = 0; j < nEntries; j ++)
	{
		chain->GetEntry(j);
		int pulser = sum->flags.pulser;
		w = sum->wais;
		double waistheta = w.theta;
		double waisphi = w.phi;
		// double theta1 = sum->mostImpulsivePeak().theta;
		// double phi1 = sum->mostImpulsivePeak().phi;
		double dist = w.distance/1e3;
		// double dtheta = theta1 - waistheta;
		// double dphi = FFTtools::wrap(phi1 - waisphi, 360, 180);
		if (treeName.EqualTo("wais")){
			if(pulser == 1){
				if(sum->deconvolved_filtered[0][0].snr<10 or sum->deconvolved_filtered[0][0].snr>14){
					continue;
				}
			// if(fabs(dtheta) > 3 or fabs(dphi)>5 ){
			// 	continue;
			// }
			// dThetavPhiH->Fill(waisphi, dtheta);
			// dThetavPhiH->Fill(sum->wais.phi, sum->peak[0][0].theta - sum->wais.theta);
			// dThetavPhiAll->Fill(sum->wais.phi, sum->peak[0][0].theta - sum->wais.theta);
			// for(int ant=0; ant < 48 ; ant++){
				// dTheta_dPhi_Power_H->Fill(FFTtools::wrap(sum->wais.phi - (sum->channels[0][ant].getPhi() - 2)*22.5,360,0), 10 - sum->wais.theta , sum->channels[0][ant].avgPower);
				// dTheta_dPhi_RMS_H->Fill(FFTtools::wrap(sum->wais.phi - (sum->channels[0][ant].getPhi() - 2)*22.5,360,0), 10 - sum->wais.theta , sum->channels[0][ant].rms);
			// }
				dThetavsdPhiH->Fill(FFTtools::wrap(sum->peak[0][0].phi - sum->wais.phi, 360, 0), sum->peak[0][0].theta - sum->wais.theta);
				dThetavsdPhi->Fill(FFTtools::wrap(sum->peak[0][0].phi - sum->wais.phi, 360, 0), sum->peak[0][0].theta - sum->wais.theta);
				dThetavsSNR->Fill(sum->deconvolved_filtered[0][0].snr, FFTtools::wrap(sum->peak[0][0].theta - sum->wais.theta, 360, 0));
				dPhivsSNR->Fill(sum->deconvolved_filtered[0][0].snr, FFTtools::wrap(sum->peak[0][0].phi - sum->wais.phi, 360, 0));
			}
			if(pulser == 5){
				if(sum->deconvolved_filtered[1][0].snr<10 or sum->deconvolved_filtered[1][0].snr>14){
					continue;
				}
				// dThetavPhiV->Fill(sum->wais.phi, sum->peak[1][0].theta - sum->wais.theta);
				// dThetavPhiAll->Fill(sum->wais.phi, sum->peak[1][0].theta - sum->wais.theta);
				// // dThetavPhiV->Fill(waisphi, dtheta);
				// for(int ant=0; ant < 48 ; ant++){
				// 	dTheta_dPhi_Power_V->Fill(FFTtools::wrap(sum->wais.phi - (sum->channels[1][ant].getPhi() - 2)*22.5,360,0), 10 - sum->wais.theta , sum->channels[1][ant].avgPower);
				// 	// dTheta_dPhi_RMS_V->Fill(FFTtools::wrap(sum->wais.phi - (sum->channels[1][ant].getPhi() - 2)*22.5,360,0), 10 - sum->wais.theta , sum->channels[1][ant].rms);
				// }
				dThetavsdPhiV->Fill(FFTtools::wrap(sum->peak[1][0].phi - sum->wais.phi, 360, 0), sum->peak[1][0].theta - sum->wais.theta);
				dThetavsdPhi->Fill(FFTtools::wrap(sum->peak[1][0].phi - sum->wais.phi, 360, 0), sum->peak[1][0].theta - sum->wais.theta);
				dThetavsSNR->Fill(sum->deconvolved_filtered[1][0].snr, FFTtools::wrap(sum->peak[1][0].theta - sum->wais.theta, 360, 0));
				dPhivsSNR->Fill(sum->deconvolved_filtered[1][0].snr, FFTtools::wrap(sum->peak[1][0].phi - sum->wais.phi, 360, 0));				
			}

		}else{
			const AnitaEventSummary::PointingHypothesis peak = sum->mostImpulsivePeak();
			dThetavsdPhi->Fill(FFTtools::wrap(peak.phi - sum->mc.phi, 360, 0), peak.theta - sum->mc.theta, sum->mc.weight);
			dThetavsSNR->Fill(peak.snr, FFTtools::wrap(peak.theta - sum->mc.theta, 360, 0), sum->mc.weight);
			dPhivsSNR->Fill(peak.snr, FFTtools::wrap(peak.phi - sum->mc.phi, 360, 0), sum->mc.weight);
		}

		// dThetavPhiAll->GetXaxis()->SetTitle("wais.phi/[degree]");
		// dThetavPhiAll->GetYaxis()->SetTitle("deltaTheta/[degree]");

		dThetavsdPhiH->GetXaxis()->SetTitle("dPhi/[degree]");
		dThetavsdPhiH->GetYaxis()->SetTitle("dTheta/[degree]");
		dThetavsdPhiV->GetXaxis()->SetTitle("dPhi/[degree]");
		dThetavsdPhiV->GetYaxis()->SetTitle("dTheta/[degree]");
		dThetavsdPhi->GetXaxis()->SetTitle("dPhi/[degree]");
		dThetavsdPhi->GetYaxis()->SetTitle("dTheta/[degree]");
		dThetavsSNR->GetXaxis()->SetTitle("snr");
		dThetavsSNR->GetYaxis()->SetTitle("dTheta/[degree]");
		dPhivsSNR->GetXaxis()->SetTitle("snr");
		dPhivsSNR->GetYaxis()->SetTitle("dPhi/[degree]");

		
	}
	TFile f2(outf.Data(), "RECREATE");
	f2.cd();
	// dThetavPhiH->Write();
	// dThetavPhiV->Write(); 
	// dThetavPhiAll->Write();
	// dTheta_dPhi_Power_H->Write();
	// dTheta_dPhi_Power_V->Write();
	// dTheta_dPhi_RMS_H->Write();
	// dTheta_dPhi_RMS_V->Write();
	dThetavsdPhiH->Write();
	dThetavsdPhiV->Write();
	dThetavsdPhi->Write();
	dThetavsSNR->Write();
	dPhivsSNR->Write();

//project 
	dThetavsdPhi->ProjectionX()->DrawNormalized()->Write();
	dThetavsdPhi->ProjectionY()->DrawNormalized()->Write();
	// dPhi->GetXaxis()->SetTitle("dPhi/[degree]");
	// dTheta->GetXaxis()->SetTitle("dTheta/[degree]");


	// dTheta_dPhi_Power_H->FitSlicesZ();
	// dTheta_dPhi_Power_V->FitSlicesZ();
	// dTheta_dPhi_RMS_H->FitSlicesZ();
	// dTheta_dPhi_RMS_V->FitSlicesZ();

	dThetavsSNR->FitSlicesY();
	dPhivsSNR->FitSlicesY();

	// gDirectory->Get("dTheta_dPhi_Power_H_1")->Write();
	// gDirectory->Get("dTheta_dPhi_Power_V_1")->Write();
	// gDirectory->Get("dTheta_dPhi_RMS_H_1")->Write();
	// gDirectory->Get("dTheta_dPhi_RMS_V_1")->Write();
	gDirectory->Get("dThetavsSNR_2")->Write();
	gDirectory->Get("dPhivsSNR_2")->Write();


	f2.Close();
}

void plotMultiSummary(){
	TString pointDir = "/Volumes/SDCard/data/wais/";
	TString fEnd = "_max_30001_sinsub_10_3_ad_2.root";	
	TString outf ="pointingSNR_Wais.root";
	int start_run = 120;
	int end_run = 155;
	TString treeName = "wais";
	_plotMultiSummary(treeName, pointDir, fEnd, outf, start_run, end_run);


	// TString pointDir = "/Volumes/SDCard/data/simulated/";
	// TString fEnd = "_max_501_sinsub_10_3_ad_2.root";	
	// TString outf ="pointingSNR_MC.root";
	// int start_run = 1;
	// int end_run = 10;
	// TString treeName = "simulation";
	// _plotMultiSummary(treeName, pointDir, fEnd, outf, start_run, end_run);

}
