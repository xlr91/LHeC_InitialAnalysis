#include "Process.h"


bool Debug = false;
bool suppress4e = false;
bool suppress4mu = false;
bool suppress2e2mu = false;

bool LepPass(GenParticle* lep_b, double ptsmear = 0){
    //Checks if the electron can be seen by the detector or not
    
    //Parameters
    double etamin = -5 ;
    double etamax = 5;
    double ptmin = 5;

    //selection criteria
    if( etamin < lep_b -> Eta && lep_b -> Eta<  etamax && ptmin < (lep_b -> PT + ptsmear)){
        return true;
    }
    return false; //lets remove the cuts
    
}  


std::vector<double> Electron_Reco(TLorentzVector scat){
    ///Electron reconstruction method
    std::vector<double> output; // Q2, x, y
    double rootS = 1300; 
    double Ee = 60; //60 GeV of initial electron
    double Qe2;
    double xe;
    double ye; 
    
    double theta = TMath::Pi() - scat.Theta();
    double Escat = abs(scat.E()); // Ee'

    ye = 1 - (Escat / Ee) * TMath::Power(  TMath::Sin(theta/2) , 2);
    Qe2 = 2 * Ee * Escat * (1 + TMath::Cos(theta));
    xe = Qe2 / (rootS * rootS * ye);
    output = {Qe2, xe, ye}; 
    
    return output;
}


std::vector<double> Hadron_Reco(TLorentzVector had){
    std::vector<double> output;
    double rootS = 1300;
    double Ee = 60;
    double Eh = abs(had.E());
    double phz = abs(had.Pz());
    double pht = abs(had.Pt());
    double Qh2;
    double xh;
    double yh; 

    yh = (Eh - phz) / (2 * Ee);
    Qh2 = pht * pht / (1-yh);
    xh = Qh2 / (rootS * rootS * yh);
    output = {Qh2, xh, yh}; 
    return output;

}


std::vector<TLorentzVector> ZZ_Reco(std::vector<GenParticle*> e_list){
    //returns a list of TLORENTZVECOTR which is [Z, Z*]
    int neg = 0;
    int pos = 2;
    int ind;

    double Zmass = 91.1876;

    std::vector<TLorentzVector> e_vecs, Ztry, ans;
    std::vector<double> Ztrymass;
    e_vecs.resize(e_list.size());
    TLorentzVector V_temp;


    for(int i = 0; i < e_list.size(); ++i){
        if(e_list[i] -> Charge == -1){
            ind = neg;
            neg++;
        }
        else if(e_list[i] -> Charge == 1){
            ind = pos;
            pos++;
        }
        V_temp.SetPtEtaPhiM(e_list[i]->PT,e_list[i]->Eta,e_list[i]->Phi,e_list[i]->Mass);
        e_vecs.at(ind) = V_temp;
    }
  
    Ztry.push_back(e_vecs[0] + e_vecs[2]);
    Ztry.push_back(e_vecs[1] + e_vecs[3]);
    Ztry.push_back(e_vecs[0] + e_vecs[3]);
    Ztry.push_back(e_vecs[1] + e_vecs[2]);

    for(int o = 0; o < Ztry.size(); ++o){
        Ztrymass.push_back(Ztry[o].M());
    }

    double curr = 0;
    int currind = 0;
    double val;

    for(int k = 0; k < Ztrymass.size(); ++k){
        val = Ztrymass[k];
        if(abs(Zmass - val) < abs(Zmass - curr)){
            curr = val;
            currind = k;
        }
    }


    switch (currind){
        case 0:
            ans.push_back(Ztry[0]);
            ans.push_back(Ztry[1]);
            break;
        case 1:
            ans.push_back(Ztry[1]);
            ans.push_back(Ztry[0]);
            break;
        case 2:
            ans.push_back(Ztry[2]);
            ans.push_back(Ztry[3]);
            break;
        case 3:
            ans.push_back(Ztry[3]);
            ans.push_back(Ztry[2]);
            break;
    }


    return ans;
}



double ptSmear(TRandom* gR, GenParticle* lep_s){
    double res;
    res = (lep_s -> PT) * (lep_s -> PT) * 8e-4;

    return gR->Gaus(0, res);
    
}

double ESmear(TRandom* gR, GenParticle* lep_s){
    double res, a, b;

    a = 15;
    b = 2;

    res = (lep_s -> E) * TMath::Sqrt( ((a * a) / (lep_s -> E)) + (b * b) ) * 1/100;
    
    return gR -> Gaus(0, res);
    
}

int main(int argc, char* argv[]) {
    // Input Delphes File

    const TString InputFile = argv[1];
    const TString OutputFileName = argv[2];

    const TString FileIdent = argv[3];
    //const TString FileIdent = 's';

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Running Process with my edit"  << std::endl;
    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "InputFile = " << InputFile << std::endl;
    std::cout << "OutputFileName = " << OutputFileName << std::endl;
    std::cout << "FileIdent = " << FileIdent << std::endl;
    std::cout << "-------------------------------------------------------------"  << std::endl;

    ExRootTreeReader * reader = NULL;
    reader = InitReader(InputFile); 

    //------------------------------------
    // Declare the output
    //------------------------------------

    OutputFile = new TFile(OutputFileName,"recreate");

    OutputFile->cd();
    OutputFile->mkdir("Example");
    OutputFile->mkdir("ParticleLevel");
    OutputFile->mkdir("4eEventLevel");
    OutputFile->mkdir("4eEventLevel/KinematicReco");
    OutputFile->mkdir("4eEventLevel/CutsAnalysis");
    OutputFile->mkdir("4eEventLevel/CutsAnalysis/Example");
    OutputFile->mkdir("4eEventLevel/Smearing");

    //Example histograms
    hEx_EventCount = new TH1D("hEx_EventCount","Event Classifications ; Event type; Number of Events",5,0,5);
    hEx_WeightCount = new TH1D("hEx_WeightCount","",10,0,1);

    hEx_Lepton_Pt = new TH1D("hEx_Lepton_Pt","Charged Lepton Events vs pT; Charged Lepton p_{T} [GeV]; Number of Particles / 5 GeV",200,0.0,1000.0);
    hEx_Z_Pt = new TH1D("hEx_Z_Pt","Z Boson Events vs pT; Z p_{T} [GeV]; Number of Particles / 5 GeV",200,0.0,1000.0);
    hEx_Jet_Pt = new TH1D("hEx_Jet_Pt","Jet Events vs pT; Jet p_{T} [GeV]; Number of Particles / 5 GeV",200,0.0,1000.0);
    hEx_ZZ_Mass = new TH1D("hEx_ZZ_Mass","Events vs mass of ZZ*; ZZ Mass [GeV]; Number of Particles / 2 GeV",125,0.0,250.0);

    //Particle-level histograms
    hPr_Jet_eta = new TH1D("hPr_Jet_eta","Jet events vs pseudorapidity; Jet #eta ; Number of Particles", 50, -6.0, 4.0);
    hPr_e_eta = new TH1D("hPr_e_eta","Electron events vs pseudorapidity; Electron #eta ; Number of Particles", 50, -7.0, 3.0);
    hPr_nue_eta = new TH1D("hPr_nue_eta","Electron Neutrino events vs pseudorapidity; Electron neutrino #eta ; Number of Particles", 50, -5.0, 5.0);
    hPr_Z_eta = new TH1D("hPr_Z_eta","Boson events vs pseudorapidity; Boson #eta ; Number of Particles", 50, -7.0, 3.0);
    hPr_H_eta = new TH1D("hPr_H_eta","Higgs events vs pseudorapidity; Higgs #eta ; Number of Particles", 50, -8.0, 2.0);

    hPr_Jet_Et = new TH1D("hPr_Jet_Et","Jet events vs transverse energy; Jet E_{T} [GeV] ; Number of Particles", 50, 0, 200);
    hPr_e_Et = new TH1D("hPr_e_Et","Electron events vs transverse energy; Electron E_{T} [GeV]; Number of Particles", 50, 0, 150);
    hPr_nue_Et = new TH1D("hPr_nue_Et","Electron Neutrino events vs transverse energy; Electron Neutrino E_{T} [GeV]; Number of Particles", 50, 0, 250);
    hPr_Z_Et = new TH1D("hPr_Z_Et","Boson vs transverse energy; Boson E_{T} [GeV] ; Number of Particles", 50, 0, 250);

    aPr_e_eta = new TEfficiency("aPr_e_eta","Acceptance of Electron vs Eta ; #eta; Acceptance",100,-10,10);
    aPr_e_Et = new TEfficiency("aPr_e_Et","Acceptance of Electron vs Transverse Energy ; E_{T} [GeV]; Acceptance", 50 , 0, 150);

    //4e event only histograms
    hEv_e_eta_nocuts = new TH1D("hEv_e_eta_nocuts","Electron (4e, no cuts) particles vs neg pseudorapidity ; Electron - #eta; Number of Particles", 50, -2.0, 8.0);
    hEv_e_eta_wicuts = new TH1D("hEv_e_eta_wicuts","Electron (4e, with cuts) particles vs pseudorapidity; Electron #eta; Number of Particles", 50, -7.0, 3.0);
    hEv_e_Et_nocuts = new TH1D("hEv_e_Et_nocuts","Electron (4e, no cuts) particles vs transverse energy; Electron E_{T}; Number of Particles", 50, 0, 150);
    hEv_e_Et_wicuts = new TH1D("hEv_e_Et_wicuts","Electron (4e, with cuts) particles vs transverse energy; Electron E_{T} ; Number of Particles", 50, 0, 150);
    hEv_nue_eta_nocuts = new TH1D("hEv_nue_eta_nocuts","Electron Neutrino (4e, no cuts) particles vs pseudorapidity; Electron neutrino #eta; Number of Particles", 50, -5.0, 5.0);
    hEv_nue_Et_nocuts = new TH1D("hEv_nue_Et_nocuts","Electron Neutrino (4e, no cuts) particles vs transverse energy; Electron Neutrino E_{T}; Number of Particles", 50, 0, 250);
    hEv_nue_eta_wicuts = new TH1D("hEv_nue_eta_wicuts","Electron Neutrino (4e, with cuts and jet check) particles vs pseudorapidity; Electron neutrino #eta; Number of Particles", 50, -5.0, 5.0);
    hEv_nue_Et_wicuts = new TH1D("hEv_nue_Et_wicuts","Electron Neutrino (4e, with cuts and jet check) particles vs transverse energy; Electron Neutrino E_{T}; Number of Particles", 50, 0, 250);
    hEv_MET_eta = new TH1D("hEv_MET_eta","Missing Energy Particle (with all cuts) vs pseudorapidity; MET #eta; Number of Particles", 50, 0, 10.0);
    hEv_MET_Et = new TH1D("hEv_MET_Et","Missing Energy Particle (with all cuts) vs transverse energy; MET E_{T}; Number of Particles", 50, 0, 250);
    aEv_e_eta = new TEfficiency("aEv_e_eta","Acceptance of Electron vs Eta (4e events); Eta; Acceptance",100,-10,10);
    aEv_e_Et = new TEfficiency("aEv_e_Et","Acceptance of Electron vs Transverse Energy (4e events); Et; Acceptance", 50 , 0, 150);
    aEv_H_eta = new TEfficiency("aEv_H_eta","Acceptance of Higgs (with jet cuts) vs Eta (4e events); Et; Acceptance", 50 , 0, 150);

    hEv_e_eta_pt = new TH2D("hEv_e_eta_pt","Electron Distribution of pt vs eta; #eta; p_{T}", 50, -7.0, 3.0, 50, 0., 150.);
    hEv_nue_eta_pt_wicuts = new TH2D("hEv_nue_eta_pt_wicuts","Electron Neutrino Distribution of pt vs eta; #eta; p_{T}", 50, -10.0, 10.0, 50, 0., 250.);
    hEv_MET_eta_pt = new TH2D("hEv_MET_eta_pt","Missing Particle Distribution of pt vs eta; #eta; p_{T}", 50, -10.0, 10.0, 50, 0., 250.);
    
    hEv_nue_MET_Phi = new TH2D("hEv_nue_MET_Phi","Phi Distribution of Missing Particle vs Electron Neutrino; #nu_{e}; Missing Particle", 50, -4.0, 4.0, 50, -4., 4.);
    hEv_nue_MET_eta = new TH2D("hEv_nue_MET_eta","Eta Distribution of Missing Particle vs Electron Neutrino; #nu_{e}; Missing Particle", 100, -5.0, 10.0, 100, -5.0, 10.0);
    hEv_nue_MET_Et = new TH2D("hEv_nue_MET_Et","Et Distribution of Missing Particle vs Electron Neutrino; #nu_{e}; Missing Particle", 50, 0, 250, 50, 0, 250);


    hEv_debugMP_Pz_E = new TH2D("hEv_debugMP_Pz_E","Checking E vs Pz distribution for Nue-MP", 100, -2000, 0, 100, 0, 2000);

    //hEv_HReco_M = new TH1D("hEv_HReco_M","Reconstructed Higgs Mass; Mass of Reconstructed Particle; Number of Particles", 50, 0, 400);
    hEv_HReco_M = new TH1D("hEv_HReco_M","Reconstructed Higgs Mass; Mass of Reconstructed Particle; Number of Particles", 100, 0, 200);
    hEv_ZReco_M = new TH1D("hEv_ZReco_M","Reconstructed Z Mass; Mass of Reconstructed Particle; Number of Particles", 100, 0, 200);
    hEv_ZstarReco_M = new TH1D("hEv_ZstarReco_M","Reconstructed Z* Mass; Mass of Reconstructed Particle; Number of Particles", 100, 0, 200);
    hEv_ZZReco_M = new TH1D("hEv_ZZReco_M","Reconstructed Z and Z* Mass; Mass of Reconstructed Particle; Number of Particles", 100, 0, 200);
    hEv_ZCurr_M = new TH1D("hEv_ZCurr_M","Reconstructed 4l Mass of Background; Mass of Reconstructed Particle (GeV); Number of Particles", 100, 0, 1000);
    
    
    hEv_jet_eta = new TH1D("hEv_jet_eta","Scattered quark  vs pseudorapidity; Jet #eta; Number of Particles", 50, -5.0, 5.0);
    hEv_H_eta = new TH1D("hEv_H_eta","Higgs vs pseudorapidity; Higgs #eta; Number of Particles", 50, -5.0, 5.0);
    
    hEvR_recoQ2_elec_hadr = new TH2D("hEvR_recoQ2_elec_hadr","2D Histogram of LogQ2, Hadron vs Electron Method; log_{10} Q2 Electron; log_{10} Q2 Hadron", 50, 0, 6.0, 50, 0., 6.);
    hEvR_recox_elec_hadr = new TH2D("hEvR_recox_elec_hadr","2D Histogram of Logx, Hadron vs Electron Method; log_{10} x Electron; log_{10} x Hadron", 50, -7, 0, 50, -7., 0.);
    hEvR_recoy_elec_hadr = new TH2D("hEvR_recoy_elec_hadr","2D Histogram of Logy, Hadron vs Electron Method; log_{10} y Electron; log_{10} y Hadron", 50, -3, 0, 50, -3., 0.);
    //hEvR_recoy_elec_hadr = new TH2D("hEvR_recoy_elec_hadr","2D Histogram of Logy, Hadron vs Electron Method; log_{10} y Electron; log_{10} y Hadron", 50, -0.3, 0.3, 50, -1.5, 0.);


    hEvR_ereco_Q2 = new TH1D("hEvR_ereco_Q2","Log Q2 values for Electron Reconstruction Method; log_{10} Q^{2}; Events", 50, 0, 6);
    hEvR_ereco_x = new TH1D("hEvR_ereco_x","Log x values for Electron Reconstruction Method; log_{10} x; Events", 50, -7, 0);
    hEvR_ereco_y = new TH1D("hEvR_ereco_y","Log y values for Electron Reconstruction Method; log_{10} y; Events", 50, -3, 0);
    hEvR_ereco_y_true = new TH1D("hEvR_ereco_y_true","y values for Electron Reconstruction Method; y; Events", 50, -1, 2);

    hEvR_hreco_Q2 = new TH1D("hEvR_hreco_Q2","Log Q2 values for Hadron Reconstruction Method; log_{10} Q^{2}; Events", 50, 0, 6);
    hEvR_hreco_x = new TH1D("hEvR_hreco_x","Log x values for Hadron Reconstruction Method; log_{10} x; Events", 50, -7, 0);
    hEvR_hreco_y = new TH1D("hEvR_hreco_y","Log y values for Hadron Reconstruction Method; log_{10} y; Events", 50, -3, 0);
    hEvR_hreco_y_true = new TH1D("hEvR_hreco_y_true"," y values for Hadron Reconstruction Method; y; Events", 50, -1, 2);
    hEvR_EPz = new TH1D("hEvR_EPz","Doing sum of E-Pz ; Sum of E-Pz; Events", 50, 0, 5000);

    hEvR_hreco_x_Q2 = new TH2D("hEvR_hreco_x_Q2","2D Histogram of Q2 vs x for hadron method; x; Q^{2}", 50, 0, 0.2, 50, 0, 45000);

    
    hEvC_Zstar = new TH1D("hEvC_Zstar","Z* Mass Cut (Less than); Z* Cut; Number of Particles below Z* cut, between 120-130", 100, 0, 100);
    hEvC_Logy = new TH1D("hEvC_Logy","Log y Cut (Less than); Log y Cut; Number of Particles below log y cut, between 120-130", 100, -3, 0);

    hEvS_e_pt = new TH1D("hEvS_e_pt","Electron particles vs transverse momentum; Electron p_{T} ; Number of Particles", 50, 0, 150);
    hEvS_e_pt_S = new TH1D("hEvS_e_pt_S","Electron smeared particles vs transverse energy; Electron p_{T} ; Number of Particles", 50, 0, 150);

    hEvS_e_E = new TH1D("hEvS_e_E","Electron particles vs Energy; Electron E ; Number of Particles", 50, 0, 500);
    hEvS_e_E_S = new TH1D("hEvS_e_E_S","Electron particles vs Energy; Electron E ; Number of Particles", 50, 0, 500);

    hEv_HReco_M_S = new TH1D("hEv_HReco_M_S","Smeared Higgs Mass; Mass of Reconstructed Particle; Number of Particles", 100, 0, 200);
    hEv_ZReco_M_S = new TH1D("hEv_ZReco_M_S","Smeared Z Mass; Mass of Reconstructed Particle; Number of Particles", 100, 0, 200);
    hEv_ZstarReco_M_S = new TH1D("hEv_ZstarReco_M_S","Smeared Z* Mass; Mass of Reconstructed Particle; Number of Particles", 100, 0, 200);

    
    


    //------------------------------------

    // Run the selection
    Process(reader, FileIdent);

    //Writes selection to the output file
    ///std::cout << "Events in EventCount: " << hEx_EventCount->GetEntries() << std::endl;
    std::cout << "Events in EventCount    : " << hEx_EventCount->GetBinContent(1) << std::endl;
    std::cout << "Events that are 4l      : " << hEx_EventCount->GetBinContent(2) << std::endl;
    std::cout << "Events seen thats 4e    : " << hEx_EventCount->GetBinContent(3) << std::endl;
    std::cout << "Events seen thats 4mu   : " << hEx_EventCount->GetBinContent(4) << std::endl;
    std::cout << "Events seen thats 2e2mu : " << hEx_EventCount->GetBinContent(5) << std::endl;
    std::cout << "Write to file..." << std::endl;

    OutputFile->cd("Example");

    TAxis *xAxis = hEx_EventCount -> GetXaxis();
    xAxis -> SetBinLabel(1, "Total Events");
    xAxis -> SetBinLabel(2, "Total 4l Events");
    xAxis -> SetBinLabel(3, "Total Observed 4e Events");
    xAxis -> SetBinLabel(4, "Total Observed 4mu Events");
    xAxis -> SetBinLabel(5, "Total Observed 2e2mu Events");

    std::cout << "Write Test 0.1" << std::endl;

    hEx_EventCount->Write();
    hEx_WeightCount->Write();

    hEx_Z_Pt->Write();
    hEx_Lepton_Pt->Write();
    hEx_Jet_Pt->Write();
    hEx_ZZ_Mass->Write();


    OutputFile->cd("ParticleLevel");

    hPr_Jet_eta -> Write();
    hPr_e_eta -> Write();
    hPr_nue_eta -> Write();
    hPr_Z_eta -> Write();
    hPr_H_eta -> Write();

    hPr_Jet_Et -> Write();
    hPr_e_Et -> Write();
    hPr_nue_Et -> Write();
    hPr_Z_Et -> Write();

    aPr_e_eta -> Write(); 
    aPr_e_Et -> Write(); 

    std::cout << "Write Test 0.2" << std::endl;


    OutputFile->cd("4eEventLevel");

    hEv_e_eta_nocuts -> Write();
    hEv_e_eta_wicuts -> Write();
    hEv_e_Et_nocuts -> Write();
    hEv_e_Et_wicuts -> Write();
    hEv_nue_eta_nocuts -> Write();
    hEv_nue_Et_nocuts -> Write();
    hEv_nue_eta_wicuts -> Write();
    hEv_nue_Et_wicuts -> Write();
    hEv_MET_eta -> Write();
    hEv_MET_Et -> Write();

    aEv_e_eta -> Write(); 
    aEv_e_Et -> Write(); 

    aEv_H_eta -> Write(); 

    hEv_e_eta_pt -> Write();
    hEv_nue_eta_pt_wicuts -> Write();
    hEv_MET_eta_pt -> Write();
    hEv_nue_MET_eta -> Write();
    hEv_nue_MET_Et -> Write();
    hEv_nue_MET_Phi -> Write();
    hEv_debugMP_Pz_E -> Write();

    

    hEv_HReco_M -> Write();
    hEv_ZReco_M -> Write();
    hEv_ZstarReco_M -> Write();
    hEv_ZZReco_M -> Write();
    hEv_ZCurr_M -> Write();


    hEv_jet_eta -> Write();
    hEv_H_eta -> Write();
  


    OutputFile->cd("4eEventLevel/KinematicReco");
    hEvR_recoQ2_elec_hadr -> Write();
    hEvR_recox_elec_hadr -> Write();
    hEvR_recoy_elec_hadr -> Write();

    hEvR_ereco_Q2 -> Write();
    hEvR_ereco_x -> Write();
    hEvR_ereco_y -> Write();
    hEvR_ereco_y_true -> Write();

    hEvR_hreco_Q2 -> Write();
    hEvR_hreco_x -> Write();
    hEvR_hreco_y -> Write();
    hEvR_hreco_y_true -> Write();

    hEvR_EPz -> Write();
    hEvR_hreco_x_Q2 -> Write();

    OutputFile->cd("4eEventLevel/CutsAnalysis");

    hEvC_Zstar -> Write();
    hEvC_Logy -> Write();

    OutputFile->cd("4eEventLevel/CutsAnalysis/Example");

    OutputFile->cd("4eEventLevel/Smearing");
    hEvS_e_pt -> Write();
    hEvS_e_pt_S -> Write();

    hEvS_e_E -> Write();
    hEvS_e_E_S -> Write();

    hEv_HReco_M -> Write();
    hEv_ZReco_M -> Write();
    hEv_ZstarReco_M -> Write();

    hEv_HReco_M_S -> Write();
    hEv_ZReco_M_S -> Write();
    hEv_ZstarReco_M_S -> Write();


	for(int n = 0; n < nCuts; ++n) {
		h_varycut.at(n)->Write();
	}

    std::cout << "Write Test 2" << std::endl;

    OutputFile->WriteObject(&cut_values, "cut_values");


    OutputFile->Close();

    std::cout << "Tidy..." << std::endl;
    delete reader;
    std::cout << "Done!" << std::endl;
    return 0;
}

ExRootTreeReader * InitReader(const TString FilePath) {

    std::cout << "InitReader" << std::endl;

    TFile * f = TFile::Open(FilePath);

    TChain * Chain = new TChain("Delphes","");

    Chain->Add(FilePath);

    // Create object of class ExRootTreeReader
    ExRootTreeReader * r = new ExRootTreeReader(Chain);
    std::cout << "InitReaderFin" << std::endl;
    return r;
}

void Process(ExRootTreeReader * treeReader, TString Ident) {

    // Get pointers to branches used in this analysis
    bEvent = treeReader->UseBranch("Event");
    bJet = treeReader->UseBranch("GenJet");
    bTruthLepton = treeReader->UseBranch("TruthLeptonParticles");
    bTruthWZ = treeReader->UseBranch("TruthWZHParticles");
    bParticle = treeReader -> UseBranch("Particle");

    Long64_t numberOfEntries = treeReader->GetEntries();
    if (Debug) numberOfEntries = 1000;

    bool is4e, is4mu, is2e2mu, good4e, good4mu, good2e2mu; //is the event a 4e event?
    bool isall4eseen, isall4museen, isall2e2museen; //can all electrons in the event be seen?
    bool goodjet; //is this jet a good jet? (Does its mass *not* correspond to a lepton)

    int nSelected = 0;
    int ecount = 0;
    int mucount = 0;

    float bscale = 8.9e-6 * 10;
    float sscale = 1.34e-5 * 10; 
    float zscale = 2.41e-6 * 100;
    float usescale;
    if(Ident == "s"){
        usescale = sscale;
    }
    else if(Ident == "b"){
        usescale = bscale;
    }
    else if(Ident == "z"){
        usescale = zscale;
    }
    else{
        usescale = 1;
    }

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Input: " << numberOfEntries << " events to process, Weight: " << usescale << std::endl;

    //Make random number generator
    gRandom = new TRandom3();
    gRandom -> SetSeed(0);
    std::vector<double> ePt_deteff, ePt_noeff, eE_deteff, eE_noeff;

    double temppts, tempEs;


    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) { 
    //for(Int_t entry = 11000; entry < 12000; ++entry) { ///CHANGE THIS BACK

        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        HepMCEvent * event = (HepMCEvent*) bEvent->At(0); 
        const float Event_Weight = usescale; 
        

        hEx_EventCount->Fill(0.5, Event_Weight);
        hEx_WeightCount->Fill(0.5,Event_Weight);

        is4e = true;
        is4mu = true;
        is2e2mu = false;
        isall4eseen = true;
        isall4museen = true;
        isall2e2museen = true;
        good4e = false;
        good4mu = false;
        good2e2mu = false;
        ecount = 0;
        mucount = 0;
        TLorentzVector Vec_Lepton_f;
        TLorentzVector Missing_Particle;
        TLorentzVector FourLepton_Vector;
        TLorentzVector FourLepton_Vector_S;
        Missing_Particle.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        FourLepton_Vector.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        FourLepton_Vector_S.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        GenParticle* my_nu;


 


        if( (entry > 0 && entry%1000 == 0) || Debug) {
            std::cout << "-------------------------------------------------------------"  << std::endl;
            std::cout << "Processing Event Number =  " << entry  << std::endl;
            std::cout << "-------------------------------------------------------------"  << std::endl;
        }

        //------------------------------------------------------------------
        // Particle Loop (Debug only)
        //------------------------------------------------------------------

        if (Debug){ // switch false with debug
            std::cout << "Event_Weight: " <<  Event_Weight << std::endl;
            for(int i = 0; i < 4; ++i){   
                GenParticle* p_Particle = (GenParticle*) bParticle->At(i);
                std::cout << "Particle " << i << " E = " << p_Particle -> E << " pZ = " << p_Particle->Pz << " pT = " << p_Particle->PT << " eta = " << p_Particle->Eta << " phi = " << p_Particle->Phi << " PID = " << p_Particle-> PID << " mass = " << p_Particle->Mass
                    << " Mother: " << p_Particle->M1 << " Daughter: " << p_Particle->D1 << " Status: " << p_Particle -> Status <<  std::endl;
            }
        }

        //------------------------------------------------------------------
        // Jet Loop
        //------------------------------------------------------------------
        for(int i = 0; i < bJet->GetEntriesFast(); ++i) {

            Jet* jet = (Jet*) bJet->At(i);

            TLorentzVector Vec_Jet;
            Vec_Jet.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);

		    hEx_Jet_Pt -> Fill( Vec_Jet.Pt(), Event_Weight);
            hPr_Jet_eta -> Fill( jet-> Eta, Event_Weight);
            hPr_Jet_Et -> Fill( TMath::Sqrt( TMath::Power(jet -> PT, 2) + TMath::Power(jet -> Mass, 2)), Event_Weight);

            if(Debug) std::cout << "Jet " << i << " pT = " << jet->PT << " eta = " << jet->Eta << " phi = " << jet->Phi << " mass = " << jet->Mass << " flavour = " << jet->Flavor << std::endl;

            //Jet Cuts
            goodjet = true;
            for(int j = 0; j < bTruthLepton->GetEntriesFast(); ++j) {
                GenParticle* lep = (GenParticle*) bTruthLepton->At(j);
                TLorentzVector Vec_Lepton_injet;
                Vec_Lepton_injet.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);
                if (Vec_Lepton_injet.DeltaR(Vec_Jet) < 0.4) goodjet = false; 
                
                if (false) std::cout << "   Jet with Lepton " << j << "  Delta R: " << Vec_Lepton_injet.DeltaR(Vec_Jet) << std::endl; //if debug
            }        
            
            if (goodjet){
                Missing_Particle = Missing_Particle - Vec_Jet;
                if (Debug){ 
                    std::cout << " MP pt: " << Missing_Particle.Pt() << " eta: " << Missing_Particle.Eta() << std::endl;
                    std::cout << " Jet E, Px, Py, Pz: " << Vec_Jet.E() << " " << Vec_Jet.Px() << " " << Vec_Jet.Py() << " " << Vec_Jet.Pz() << std::endl;
                    std::cout << " MP  E, Px, Py, Pz: " << Missing_Particle.E() << " " << Missing_Particle.Px() << " " << Missing_Particle.Py() << " " << Missing_Particle.Pz() << std::endl; 
                }
            }

        } // Jet Loop
        

        //------------------------------------------------------------------
        // Lepton Loops
        //------------------------------------------------------------------

        //first pass to check if the particles/event is good or not
        for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) {

            GenParticle* lep = (GenParticle*) bTruthLepton->At(i);

            if(Debug) {
                std::cout << "Lepton " << i << " pT = " << lep->PT << " eta = " << lep->Eta << " phi = " << lep->Phi << " PID = " << lep-> PID << " mass = " << lep->Mass
                << " Mother: " << lep->M1 <<  std::endl;
            }
  
            TLorentzVector Vec_Lepton;
            Vec_Lepton.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);
            

            // Look for electrons or muons
            if( abs(lep->PID) == 11 || abs(lep->PID) == 13 ) {
                hEx_Lepton_Pt->Fill( Vec_Lepton.Pt(), Event_Weight );  
            }

            if( abs(lep->PID) == 11) {
                temppts = ptSmear(gRandom, lep);
                tempEs = ESmear(gRandom, lep);

                ePt_deteff.push_back(temppts);
                ePt_noeff.push_back(lep -> PT);
                
                eE_deteff.push_back(tempEs);
                eE_noeff.push_back(lep->E);
                
                if(Debug) std::cout<<"smearing pt: "<< temppts << std::endl;
                hPr_e_eta -> Fill( lep-> Eta ,Event_Weight);
                hPr_e_Et -> Fill( TMath::Sqrt( TMath::Power(lep -> PT, 2) + TMath::Power(lep -> Mass, 2)), Event_Weight);
                
                aPr_e_eta -> FillWeighted(LepPass(lep, temppts), Event_Weight, lep->Eta);
                aPr_e_Et -> FillWeighted(LepPass(lep, temppts), Event_Weight, TMath::Sqrt( TMath::Power(lep -> PT, 2) + TMath::Power(lep -> Mass, 2)));

                //check if all 4e events seen 
                if(!LepPass(lep, temppts)){ isall4eseen = false; isall2e2museen = false;}
                is4mu = false;  
                ecount++;


                
            }

            if( abs(lep->PID) == 12) {
                if (Debug) std::cout << " Nue E, Px, Py, Pz, P: " << Vec_Lepton.E() << " " << Vec_Lepton.Px() << " " << Vec_Lepton.Py() << " " << Vec_Lepton.Pz() << " " << Vec_Lepton.P()<<  std::endl;
                if (Debug) std::cout << " Nue theta, sintheta, costheta: "<<Vec_Lepton.Theta() <<  " " << TMath::Sin(Vec_Lepton.Theta()) << " " << TMath::Cos(Vec_Lepton.Theta()) <<std::endl;
                hPr_nue_eta -> Fill( lep-> Eta ,Event_Weight);
                hPr_nue_Et -> Fill( TMath::Sqrt( TMath::Power(lep -> PT, 2) + TMath::Power(lep -> Mass, 2)), Event_Weight);
            }

            //if it contains a muon, its definitely not a 4e event 
            if (abs(lep -> PID) == 13){
                if(!LepPass(lep)) {isall4museen = false ;  isall2e2museen = false;}  
                is4e = false;
                mucount++;

                temppts = ptSmear(gRandom, lep);
                tempEs = ESmear(gRandom, lep);

                ePt_deteff.push_back(temppts);
                ePt_noeff.push_back(lep -> PT);
                eE_deteff.push_back(tempEs);
                eE_noeff.push_back(lep->E);
                
            }
        } // Lepton Loop

        if(ecount == 2 && mucount == 2) is2e2mu = true;
        
        if(suppress4e) is4e = false;
        if(suppress4mu) is4mu = false;
        if(suppress2e2mu) is2e2mu = false;

        

        if (is4e || is4mu || is2e2mu){
            TLorentzVector Debug_MP;
            std::vector<GenParticle*> e_par_list;
            if (Debug) {
                if(is4e) std::cout << "This is a 4e event. Seen: " << isall4eseen << std::endl;  
                if(is4mu) std::cout << "This is a 4mu event. Seen: " << isall4museen << std::endl;  
                if(is2e2mu) std::cout << "This is a 2e2mu event. Seen: " << isall2e2museen << std::endl;  
            }
            

            if ((isall4eseen || isall4museen || isall2e2museen) && Missing_Particle.Eta() != 0) hEv_jet_eta -> Fill(Missing_Particle.Eta(), Event_Weight);
            //into lepton for loop again, this time its only  4e events
            bool interestingjet = true;
            if(Missing_Particle.Eta() == 0) interestingjet = false;


            int lepint_temp = 0;
            for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i){
                GenParticle* lep_e = (GenParticle*) bTruthLepton->At(i);
                TLorentzVector Vec_Lepton_e;
                TLorentzVector Vec_Lepton_e_S;

                //std::cout<<lepint_temp <<std::endl;
                Vec_Lepton_e.SetPtEtaPhiM(lep_e->PT,lep_e->Eta,lep_e->Phi,lep_e->Mass); //change this for smearing
                

                if( abs(lep_e->PID) == 12) {
                    hEv_nue_eta_nocuts -> Fill( lep_e-> Eta ,Event_Weight);
                    hEv_nue_Et_nocuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);
                    if (isall4eseen || isall4museen || isall2e2museen){
                        my_nu = lep_e;
                        //fill the electron neutrino comparison with missing energy
                        if(interestingjet) hEv_nue_eta_wicuts -> Fill( -(lep_e-> Eta) ,Event_Weight);
                        hEv_nue_Et_wicuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);

                        hEv_nue_eta_pt_wicuts -> Fill(lep_e-> Eta,  lep_e-> PT);
                    }
                }
                
                if( abs(lep_e->PID) == 11 || abs(lep_e->PID) == 13){
                    hEv_e_eta_nocuts -> Fill( -(lep_e -> Eta) , Event_Weight);
                    hEv_e_Et_nocuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);

                    //hEvS_e_pt_S -> Fill(TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)) + ePt_deteff.at(lepint_temp), Event_Weight);
                    
                    
                    
                    

                    if(Debug){
                        std::cout << "Lepton Number " << i << "vector loc" << lepint_temp <<  std::endl;
                        std::cout << "smeared pT: " << lep_e -> PT  + ePt_deteff.at(lepint_temp) << std::endl;
                        std::cout << "smeared E: " << lep_e -> E << " " << eE_deteff.at(lepint_temp) << std::endl;
                    } 
                

                    aEv_e_eta -> FillWeighted(LepPass(lep_e, ePt_deteff.at(lepint_temp)), Event_Weight, lep_e->Eta);
                    aEv_e_Et -> FillWeighted(LepPass(lep_e, ePt_deteff.at(lepint_temp)), Event_Weight, TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)));

                    if(LepPass(lep_e, ePt_deteff.at(lepint_temp))){
                        hEv_e_eta_wicuts -> Fill( lep_e -> Eta , Event_Weight);
                        hEv_e_Et_wicuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);   
                        
                        e_par_list.push_back(lep_e);

                        //Vec_Lepton_e_S.SetPtEtaPhiE(lep_e -> PT  + ePt_deteff.at(lepint_temp) ,lep_e->Eta,lep_e->Phi, lep_e->E + eE_deteff.at(lepint_temp)); //change this for smearing (ORIGINAL)
                        

                        
                        if(abs(lep_e->PID) == 11){
                            Vec_Lepton_e_S.SetPxPyPzE(lep_e -> Px + eE_deteff.at(lepint_temp), lep_e -> Py + eE_deteff.at(lepint_temp), lep_e -> Pz + eE_deteff.at(lepint_temp), lep_e -> E + eE_deteff.at(lepint_temp)); //change this for smearing (NEW)
                        }

                        if(abs(lep_e->PID) == 13){
                            Vec_Lepton_e_S.SetPxPyPzE(lep_e -> Px + ePt_deteff.at(lepint_temp), lep_e -> Py + ePt_deteff.at(lepint_temp), lep_e -> Pz + ePt_deteff.at(lepint_temp), lep_e -> E + ePt_deteff.at(lepint_temp)); //change this for smearing (NEW)
                        }

                        
                        Missing_Particle = Missing_Particle - Vec_Lepton_e;
                        FourLepton_Vector = FourLepton_Vector + Vec_Lepton_e;
                        FourLepton_Vector_S = FourLepton_Vector_S + Vec_Lepton_e_S;
                    }

                    hEv_e_eta_pt -> Fill(lep_e -> Eta, lep_e -> PT);
                    lepint_temp++;
                }

            }

            hEx_EventCount -> Fill(1.5, Event_Weight);

            
            // Event santiy check here. It needs to be 4 events, and it needs to be 2 pos and 2 neg
            if (e_par_list.size() != 4){
                isall4eseen = false;
                isall4museen = false;
                isall2e2museen = false;
            } 

            int sanityp = 0;
            int sanityn = 0;

            for(int i = 0; i < e_par_list.size(); ++i){
                if(e_par_list[i] -> Charge == -1){
                    sanityn++;
                }
                else if(e_par_list[i] -> Charge == 1){
                    sanityp++;
                }
            }

            if(!(sanityp == 2 && sanityn == 2)){
                isall4eseen = false;
                isall4museen = false;
                isall2e2museen = false;
            }
            
            if(isall4eseen && is4e) good4e = true;
            if(isall4museen && is4mu) good4mu = true;
            if(isall2e2museen && is2e2mu) good2e2mu = true;
            if(Debug) std::cout << good4e << good4mu << std::endl;


            if (good4e || good4mu || good2e2mu) {
                if(is4e)  hEx_EventCount -> Fill(2.5, Event_Weight);
                if(is4mu) hEx_EventCount -> Fill(3.5, Event_Weight);
                if(is2e2mu) hEx_EventCount -> Fill(4.5, Event_Weight);
                
                TLorentzVector my_nu_vec;
                
                my_nu_vec.SetPtEtaPhiM(my_nu->PT,my_nu->Eta,my_nu->Phi,my_nu->Mass);
                
                Debug_MP =   my_nu_vec - Missing_Particle;

                if(is4e&&Debug) std::cout << "4e" << std::endl;
                if(is4mu&&Debug) std::cout << "4mu" << std::endl;

                TLorentzVector Hadronic_Vector = -Missing_Particle;

                std::vector<double> E_Reco_Lst = Electron_Reco(my_nu_vec);
                std::vector<double> H_Reco_Lst = Hadron_Reco(Hadronic_Vector);
                std::vector<TLorentzVector> Z_Reco_Lst = ZZ_Reco(e_par_list);  

                double sumEPz = my_nu_vec.E() - my_nu_vec.Pz() + ( Hadronic_Vector.E() - Hadronic_Vector.Pz());
                
                //Z* and logycut cut analysis// My Version
                for(int zzz = 1; zzz <= hEvC_Zstar -> GetNbinsX(); zzz++){
                    if( 120 <= FourLepton_Vector.M() && FourLepton_Vector.M() <= 130){
                        if(Z_Reco_Lst[1].M() < hEvC_Zstar->GetBinLowEdge(zzz)){
                            hEvC_Zstar -> Fill(hEvC_Zstar -> GetBinCenter(zzz), Event_Weight);
                        }
                        if(TMath::Log10(H_Reco_Lst[2]) < hEvC_Logy -> GetBinLowEdge(zzz)){
                            hEvC_Logy -> Fill(hEvC_Logy -> GetBinCenter(zzz), Event_Weight);
                        }
                    }
                }

                if (Debug){
                    std::cout << "MET pt: " << Missing_Particle.Pt()  << " eta: " << Missing_Particle.Eta() << std::endl;
                    std::cout << "E-Reco Q2: " << E_Reco_Lst[0] << " x: " << E_Reco_Lst[1] << " y: " << E_Reco_Lst[2] << std::endl;
                    std::cout << "H-Reco Q2: " << H_Reco_Lst[0] << " x: " << H_Reco_Lst[1] << " y: " << H_Reco_Lst[2] << std::endl;
                    std::cout << "Reco Higgs Mass?: " << FourLepton_Vector.M() << std::endl;
                    std::cout << "SumPz: " << sumEPz << std::endl;

                    std::cout << "ELIST CHECK " << e_par_list[0] -> Charge  << std::endl;
                    
                    for(int oo = 0; oo < 2; ++oo){
                        std::cout << Z_Reco_Lst[oo].M() << std::endl;
                    }

                    //std::cout << "Random Number: " << ptSmear(gRandom) << std::endl;

                }



                
                hEv_MET_eta -> Fill(Missing_Particle.Eta() ,Event_Weight);
                hEv_MET_Et -> Fill(Missing_Particle.Pt() ,Event_Weight); //pt = et for neutrino
                
                hEv_MET_eta_pt -> Fill(Missing_Particle.Eta(), Missing_Particle.Pt());
                hEv_nue_MET_eta -> Fill(my_nu->Eta, Missing_Particle.Eta());
                hEv_nue_MET_Et -> Fill(my_nu -> PT, Missing_Particle.Pt());
                hEv_nue_MET_Phi -> Fill(my_nu -> Phi, Missing_Particle.Phi());

                hEv_debugMP_Pz_E -> Fill(Debug_MP.Pz(), Debug_MP.E());

                hEv_HReco_M -> Fill(FourLepton_Vector.M(), Event_Weight); //not the hadronic vector 
                hEv_HReco_M_S -> Fill(FourLepton_Vector_S.M(), Event_Weight); //not the hadronic vector 
                hEv_ZReco_M -> Fill(Z_Reco_Lst[0].M(), Event_Weight);
                hEv_ZCurr_M -> Fill(FourLepton_Vector.M(), Event_Weight);
                hEv_ZstarReco_M -> Fill(Z_Reco_Lst[1].M(), Event_Weight);

                hEv_ZZReco_M -> Fill(Z_Reco_Lst[0].M(), Event_Weight);
                hEv_ZZReco_M -> Fill(Z_Reco_Lst[1].M(), Event_Weight);

                hEv_H_eta -> Fill(-(FourLepton_Vector.Eta()), Event_Weight);
    


                //try{
                hEvR_recoQ2_elec_hadr -> Fill(TMath::Log10(E_Reco_Lst[0]), TMath::Log10(H_Reco_Lst[0]));
                hEvR_recox_elec_hadr -> Fill(TMath::Log10(E_Reco_Lst[1]), TMath::Log10(H_Reco_Lst[1]));
                hEvR_recoy_elec_hadr -> Fill(TMath::Log10(E_Reco_Lst[2]), TMath::Log10(H_Reco_Lst[2]));
                
                hEvR_ereco_Q2 -> Fill(TMath::Log10(E_Reco_Lst[0]), Event_Weight);
                hEvR_ereco_x -> Fill(TMath::Log10(E_Reco_Lst[1]), Event_Weight);
                hEvR_ereco_y -> Fill(TMath::Log10(E_Reco_Lst[2]), Event_Weight);
                hEvR_ereco_y_true -> Fill(E_Reco_Lst[2], Event_Weight);
                
                hEvR_hreco_Q2 -> Fill(TMath::Log10(H_Reco_Lst[0]), Event_Weight);
                hEvR_hreco_x -> Fill(TMath::Log10(H_Reco_Lst[1]), Event_Weight);
                hEvR_hreco_y -> Fill(TMath::Log10(H_Reco_Lst[2]), Event_Weight);
                hEvR_hreco_y_true -> Fill(H_Reco_Lst[2], Event_Weight);
                hEvR_hreco_x_Q2 -> Fill(H_Reco_Lst[1], H_Reco_Lst[0]);

                hEvR_EPz -> Fill(sumEPz, Event_Weight);


                
                for(int indx=0; indx < 4; ++indx){
                    hEvS_e_pt -> Fill(ePt_noeff.at(indx), Event_Weight);
                    hEvS_e_pt_S -> Fill(ePt_noeff.at(indx) + ePt_deteff.at(indx), Event_Weight); 
                    hEvS_e_E -> Fill(eE_noeff.at(indx), Event_Weight);
                    hEvS_e_E_S ->Fill(eE_noeff.at(indx) + eE_deteff.at(indx), Event_Weight); 
                    if(Debug) std::cout << eE_deteff.at(indx) << std::endl;
                }
                
                    
            }
            
        }
	            
        //------------------------------------------------------------------
        // W/Z/H Boson Loop
        //------------------------------------------------------------------
        std::vector<TLorentzVector> list_Zboson;
        for(int i = 0; i < bTruthWZ->GetEntriesFast(); ++i) {
            
            GenParticle* v = (GenParticle*) bTruthWZ->At(i);
            
            if(Debug) std::cout << "Boson " << i << " pT = " << v->PT << " eta = " << v->Eta << " phi = " << v->Phi << " PID = " << v-> PID << " mass = " << v->Mass 
                << " Daughter: " << v->D1 <<  std::endl;
        
        
            TLorentzVector Vec_V;
            Vec_V.SetPtEtaPhiM(v->PT,v->Eta,v->Phi,v->Mass);


            
            if (v -> Mass == 125){
                hPr_H_eta -> Fill( v->Eta, Event_Weight);
            }
            

    		// Look for Z bosons (ID = 23)
            if( v->PID == 23 ) {
                hEx_Z_Pt->Fill( Vec_V.Pt(), Event_Weight );
                list_Zboson.push_back(Vec_V);
                hPr_Z_eta->Fill( v->Eta, Event_Weight );
                hPr_Z_Et -> Fill( TMath::Sqrt( TMath::Power(v -> PT, 2) + TMath::Power(v -> Mass, 2)), Event_Weight);
                
            }
            
        } // W/Z/H Loop

	    if( list_Zboson.size() > 1 ) {

		    TLorentzVector Higgs = list_Zboson.at(0) + list_Zboson.at(1);

            hEx_ZZ_Mass->Fill( Higgs.M() , Event_Weight );
            aEv_H_eta ->  FillWeighted(isall4eseen && is4e, Event_Weight, Higgs.Eta());
            
            //do the higgs selectoin events here 
            //as in acceptance of the higgs eta and pt
	    }

        ePt_deteff.clear();
        ePt_noeff.clear();
        eE_deteff.clear();
        eE_noeff.clear();
    } // Loop over all events
}

