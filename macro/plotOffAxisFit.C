#include "TFile.h"
#include "TTree.h"
#include "AnitaEventSummary.h" 
#include "AllCorrelationSummaryAnita4.h"
#include "AntennaPositions.h" 

#include "FFTtools.h"
// plot delta theta vs phi for all wais event runs.
void plotOffAxisFit()
{
	TString PlotPrefix = "Plot_";
	TString pointDir = "../drivers/corrTrees/run";
	// TString fEnd = "_max_30002_sinsub_10_3_ad_2.root";	
	TString fEnd = "CorrTree_4.root";	
	TString outf = pointDir + PlotPrefix + fEnd;
	const UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance(); 
	
	TChain* chain = new TChain("corrTree");
	for(int i = 123; i < 153; i++)
	// for(int i = 139; i < 141; i++)
	{
		if(i==129 || i==130 || i==131 || i==132) continue;
		//if(i==132) continue;
		TString inf = pointDir + TString::Itoa(i, 10) + fEnd;
		chain->Add(inf);
	}


	TFile *file = TFile::Open("../drivers/corrTrees/allRuns" + fEnd,"RECREATE");
	chain->CloneTree(-1,"fast");
	file->Write();

	Int_t pol = 0;
	AllCorrelationSummaryAnita4* corr = new AllCorrelationSummaryAnita4;
	chain->SetBranchAddress("corr", &corr);
	// chain->SetBranchAddress("pol", &pol);

	Long_t nEntries = chain->GetEntries();

	TH2D* GroupDelayvsPhiH = new TH2D("GroupDelayvsPhiH", "GroupDelayvsPhiH", 100, -45,45, 100, -0.26, 0.26);	
	TH2D* GroupDelayvsPhiV = new TH2D("GroupDelayvsPhiV", "GroupDelayvsPhiV", 100, -45,45, 100, -0.26, 0.26);	
	TH2D* GroupDelayvsPhiAll = new TH2D("GroupDelayvsPhiAll", "GroupDelayvsPhiAll", 100, -45,45, 100, -0.26, 0.26);	

	GroupDelayvsPhiH->GetXaxis()->SetTitle("PhiWave - MiddlePhi/[degree]");
	GroupDelayvsPhiH->GetYaxis()->SetTitle("#Delta#Deltat/[ns]");
	GroupDelayvsPhiV->GetXaxis()->SetTitle("PhiWave - MiddlePhi/[degree]");
	GroupDelayvsPhiV->GetYaxis()->SetTitle("#Delta#Deltat/[ns]");
	GroupDelayvsPhiAll->GetXaxis()->SetTitle("PhiWave - MiddlePhi/[degree]");
	GroupDelayvsPhiAll->GetYaxis()->SetTitle("#Delta#Deltat/[ns]");


	for (int j = 0; j < nEntries; j ++)
	{
		chain->GetEntry(j);
		pol = corr->pol;


		 int pair = 0;
	   //Now we can make correlations an
		for(int i = 0; i < 14; i++){
			for(int j = i + 1; j < 15; j++){
				int ant1 = corr->fifteenAnts[i];
				int ant2 = corr->fifteenAnts[j];
				int sign = 1;
				if((ant2-ant1+16+2)%16>4){  // do not consider ant pair that separate larger than 2 phi sectors
				  continue;
				} 
				// if(1){// all pairs
				// if((ant2-ant1+16+2)%16==2){ // same phi sector
				if((ant2-ant1+16+2)%16==1  or (ant2-ant1+16+2)%16==3){  // 1 phi sector
				// if((ant2-ant1+16+2)%16==0  or (ant2-ant1+16+2)%16==4){  // 2 phi sector
					if((ant2-ant1+16+2)%16==0 or (ant2-ant1+16+2)%16==1) {
						sign = -1;
					}else{
						sign = 1;
					}
					double measuredDeltaT = corr->maxCorTimes[pair];
					GroupDelayvsPhiAll->Fill(FFTtools::wrap(corr->phiWave*180/3.14 - (ap->phiAnt[pol][ant1] + ap->phiAnt[pol][ant2])/2.0, 180 , 0) , sign*(measuredDeltaT - corr->expectedDeltaT[pair]));
					if(pol == 0){
						GroupDelayvsPhiH->Fill(FFTtools::wrap(corr->phiWave*180/3.14 - (ap->phiAnt[pol][ant1] + ap->phiAnt[pol][ant2])/2.0, 180 , 0) , sign*(measuredDeltaT - corr->expectedDeltaT[pair]));
					}
					if(pol == 1){
						GroupDelayvsPhiV->Fill(FFTtools::wrap(corr->phiWave*180/3.14 - (ap->phiAnt[pol][ant1] + ap->phiAnt[pol][ant2])/2.0, 180 , 0) , sign*(measuredDeltaT - corr->expectedDeltaT[pair]));
					}
				}
			 	pair++;
			}
		}

		// std::cout<< pol<<" "<< std::endl;

		// for(int pair=0; pair < 49 ; pair++){

			
		// 	int ant1 = corr->firstAnt[pair];
		// 	int ant2 = corr->secondAnt[pair];
		// 	if(ant1==45 or ant2 == 45) continue;

		// 	// if((pair<=5) or (pair>=36 and pair<=40)){   // only  same phi sector 
		// 	// if((pair>=6 and pair<=35) or (pair>=41)){  // only different phi sector 
		// 	// if((pair>=6 and pair<=19) or(pair>=26 and pair<=35) or (pair>=41)){  // only 1 phi sector 
		// 	if(pair>=20 and pair<=25){  // only 2 phi sector 
		// 		double measuredDeltaT = corr->maxCorTimes[pair];
		// 		GroupDelayvsPhiAll->Fill(FFTtools::wrap(corr->phiWave*180/3.14 - (ap->phiAnt[pol][ant1] + ap->phiAnt[pol][ant2])/2.0, 180 , 0) , measuredDeltaT - corr->expectedDeltaT[pair]);
		// 		if(pol == 0){
		// 			GroupDelayvsPhiH->Fill(FFTtools::wrap(corr->phiWave*180/3.14 - (ap->phiAnt[pol][ant1] + ap->phiAnt[pol][ant2])/2.0, 180 , 0) , measuredDeltaT - corr->expectedDeltaT[pair]);
		// 		}
		// 		if(pol == 1){
		// 			GroupDelayvsPhiV->Fill(FFTtools::wrap(corr->phiWave*180/3.14 - (ap->phiAnt[pol][ant1] + ap->phiAnt[pol][ant2])/2.0, 180 , 0) , measuredDeltaT - corr->expectedDeltaT[pair]);
		// 		}
		// 	}
		// }
	}

	

	TFile f2(outf.Data(), "RECREATE");
	f2.cd();

	GroupDelayvsPhiH->Write();
	GroupDelayvsPhiV->Write();
	GroupDelayvsPhiAll->Write();
	GroupDelayvsPhiH->FitSlicesY();
	GroupDelayvsPhiV->FitSlicesY();
	GroupDelayvsPhiAll->FitSlicesY();
	gDirectory->Get("GroupDelayvsPhiH_1")->Write();
	gDirectory->Get("GroupDelayvsPhiV_1")->Write();
	gDirectory->Get("GroupDelayvsPhiAll_1")->Write();



	// GroupDelayvsPhiH->Write();
	// GroupDelayvsPhiV->Write();
	// GroupDelayvsPhiAll->Write();

	// GroupDelayvsPhiH->FitSlicesY();
	// GroupDelayvsPhiV->FitSlicesY();
	// GroupDelayvsPhiAll->FitSlicesY();
	// gDirectory->Get("GroupDelayvsPhiH_1")->Write();
	// gDirectory->Get("GroupDelayvsPhiV_1")->Write();
	// gDirectory->Get("GroupDelayvsPhiAll_1")->Write();


	// dTheta_dPhi_Power_H->FitSlicesZ();
	// dTheta_dPhi_Power_V->FitSlicesZ();
	// gDirectory->Get("dTheta_dPhi_Power_H_1")->Write();
	// gDirectory->Get("dTheta_dPhi_Power_V_1")->Write();


	f2.Close();
}

