#include "Process.h"




bool Debug = false;

bool LepPass(GenParticle* lep_b){
    //Checks if the electron can be seen by the detector or not
    
    //Parameters
    double etamin = -5 ;
    double etamax = 5;
    double ptmin = 5;

    //selection criteria
    if( etamin < lep_b -> Eta && lep_b -> Eta<  etamax && ptmin < lep_b -> PT){
        return true;
    }
    return false; //lets remove the cuts
    //return true; //lets remove the cuts
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

    std::vector<GenParticle*> e_sorted;
    std::vector<TLorentzVector> e_vecs, Ztry, ans;
    std::vector<double> Ztrymass;
    e_sorted.resize(e_list.size());
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
        e_sorted[ind] = e_list[i];
        V_temp.SetPtEtaPhiM(e_list[i]->PT,e_list[i]->Eta,e_list[i]->Phi,e_list[i]->Mass);
        e_vecs[ind] = V_temp;
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

int main(int argc, char* argv[]) {

    // Input Delphes File

    const TString InputFile = argv[1];
    const TString OutputFileName = argv[2];
    const TString FileIdent = argv[3];

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

    //Example histograms
    hEx_EventCount = new TH1D("hEx_EventCount","Event Classifications ; Event type; Number of Events",4,0,4);
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
    hEv_HReco_M = new TH1D("hEv_HReco_M","Reconstructed Higgs Mass; Mass of Reconstructed Particle; Number of Particles", 75, 0, 150);
    hEv_ZReco_M = new TH1D("hEv_ZReco_M","Reconstructed Z Mass; Mass of Reconstructed Particle; Number of Particles", 75, 0, 150);
    hEv_ZstarReco_M = new TH1D("hEv_ZstarReco_M","Reconstructed Z* Mass; Mass of Reconstructed Particle; Number of Particles", 75, 0, 150);
    hEv_ZZReco_M = new TH1D("hEv_ZZReco_M","Reconstructed Z and Z* Mass; Mass of Reconstructed Particle; Number of Particles", 75, 0, 150);
    
    
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
    
    
    //logyCut // Andy's Example // Initializations
    
    nCuts = 100;
    min_cut = -3;
    double step_cut =  abs(min_cut) / nCuts;
    

    for(int n=0; n<nCuts; ++n){
        TString suffix = "_cut";
        suffix += n;
        TString titlething = Form("%d",min_cut + n*step_cut);

        TH1D* h_new_cuts = new TH1D("h_varycut"+suffix, titlething + "; m_{4l} [GeV]; Events / 4 GeV", 50, 50, 250);
        h_varycut.push_back(h_new_cuts);
        cut_values.push_back(min_cut + n*step_cut);
        if(Debug) std::cout<<min_cut + n*step_cut << std::endl;
    }

    //------------------------------------

    // Run the selection
    Process(reader, FileIdent);

    //Writes selection to the output file
    ///std::cout << "Events in EventCount: " << hEx_EventCount->GetEntries() << std::endl;
    std::cout << "Events in EventCount: " << hEx_EventCount->GetBinContent(1) << std::endl;
    std::cout << "Events that are 4e  : " << hEx_EventCount->GetBinContent(2) << std::endl;
    std::cout << "Events seen thats 4e: " << hEx_EventCount->GetBinContent(3) << std::endl;
    std::cout << "Events seen 4e GJet : " << hEx_EventCount->GetBinContent(4) << std::endl;
    std::cout << "Write to file..." << std::endl;

    OutputFile->cd("Example");

    TAxis *xAxis = hEx_EventCount -> GetXaxis();
    xAxis -> SetBinLabel(1, "Total Events");
    xAxis -> SetBinLabel(2, "Total 4e Events");
    xAxis -> SetBinLabel(3, "Total 4e Events Passing Cuts");
    xAxis -> SetBinLabel(4, "Total 4e passed cuts with Good Jets");

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

    //hEv_HReco_M->SetOption("L");
    //hEv_HReco_M->SetFillColorAlpha (2, 1);
    //hEv_HReco_M->SetLineColor(3);
    //std::cout<<"Yeet" << std::endl;
    //std::cout<<hEv_HReco_M -> GetOption() << std::endl;
    //hEv_HReco_M -> GetXaxis()->SetRangeUser(100, 150);
    hEv_HReco_M -> Write();
    hEv_ZReco_M -> Write();
    hEv_ZstarReco_M -> Write();
    hEv_ZZReco_M -> Write();
    /*
        TCanvas *c1 = new TCanvas("c1"," Test",50,50,1680,1050); //initialize TCanvas object
        TLegend *legend1;  //For now, just defines the legend
        double maxhist = hEv_HReco_M -> GetMaximum();
        if(maxhist < hEv_ZReco_M -> GetMaximum()) maxhist = hEv_ZReco_M -> GetMaximum();
        if(maxhist < hEv_ZstarReco_M -> GetMaximum()) maxhist = hEv_ZstarReco_M -> GetMaximum();

        
        //note that here we're applying the formatting on the _first_ histogram we wanna stick on tcanvas
        hEv_HReco_M -> SetAxisRange(0, maxhist*1.1, "Y");
        hEv_HReco_M -> SetTitle("Reconstruction of Higgs, Z, and Z* Bosons"); //sets title of the histogram ur gonna print
        hEv_HReco_M ->GetXaxis()->SetTitle("Mass (GeV)");  //sets x axis title
        hEv_HReco_M ->GetYaxis()->SetTitle("Number of Events"); //sets y axis title
        hEv_HReco_M -> GetXaxis()->SetRangeUser(0, 130); //sets the visual range, like my old histogram goes from 0 to 150 but i wanna show till 130
        hEv_HReco_M -> SetLineColor(kRed); //changes the bars of the histogram so its red
        hEv_HReco_M -> SetStats(kFALSE); //removes the annoying 'mean median and entries' box thing
        hEv_HReco_M -> Draw("hist E2"); //Draws the histogram using the option hist and E2
        //The Draw command takes hEv_HReco_M with its applied options and sticks it into the TCanvas canvas

        hEv_ZReco_M -> SetLineColor(kBlue); //this is my second histogram
        hEv_ZReco_M -> Draw("hist same E2"); //NOTE THE 'SAME' BIT, draws the histogram on the same canvas

        hEv_ZstarReco_M->SetLineColor(kGreen); //third histogram
        hEv_ZstarReco_M -> Draw("same hist E2");

        legend1 = new TLegend(0.1,0.8,0.25,0.9); // Initialize the TLegend object that you place within the TCanvas object
        //numbers are (x1, y1, x2, y2), bottom left and top right corner (i think)
        legend1->SetHeader("Particles","C"); // option "C" allows to center the header. 'Particles' is the Title of the legend
        legend1->AddEntry(hEv_HReco_M, "Higgs"); //Adds the entry and the label
        legend1->AddEntry(hEv_ZReco_M,"Z Boson");
        legend1->AddEntry(hEv_ZstarReco_M,"Z* Boson");
        legend1-> Draw("same");  //Draws on same histogram

        c1 -> Write("ReconstructedMass"); //This final command saves the canvas that we've made art on into ur root file
        //c1 -> Print("TestCanvas.png"); //This command prints the TCanvas as a png which you can view immediately on vscode, highly suggested
    */


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

	for(int n = 0; n < nCuts; ++n) {
		h_varycut.at(n)->Write();
	}

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

    bool is4e; //is the event a 4e event?
    bool isalleseen; //can all electrons in the event be seen?
    bool goodjet; //is this jet a good jet? (Does its mass *not* correspond to a lepton)

    int nSelected = 0;

    float bscale = 8.9e-6 * 10;
    float sscale = 1.34e-5 * 10; 
    float usescale;
    if(Ident == "s"){
        usescale = sscale;
    }
    else if(Ident == "b"){
        usescale = bscale;
    }

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Input: " << numberOfEntries << " events to process, Weight: " << usescale << std::endl;






    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) { 
    //for(Int_t entry = 11000; entry < 12000; ++entry) { ///CHANGE THIS BACK


        //if (entry == 1000) Debug = false;

        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        HepMCEvent * event = (HepMCEvent*) bEvent->At(0); 
    	//const float Event_Weight = event->Weight;
        const float Event_Weight = usescale; //THIS IS WHERE THE EVENT WEIGHT NEEDS TO BE CHANGED FOR THE BACKGROUND SAMPLE AND THE SIGNAL
        //const float Event_Weight = 1;  //TESTING PURPOSE, REMOVE THIS WHEN DONE

        hEx_EventCount->Fill(0.5, Event_Weight);
        hEx_WeightCount->Fill(0.5,Event_Weight);

        is4e = true;
        isalleseen = true;
        TLorentzVector Vec_Lepton_f;
        TLorentzVector Missing_Particle;
        TLorentzVector FourLepton_Vector;
        Missing_Particle.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        FourLepton_Vector.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        GenParticle* my_nu;

        if( (entry > 0 && entry%1000 == 0) || Debug) {
            std::cout << "-------------------------------------------------------------"  << std::endl;
            std::cout << "Processing Event Number =  " << entry  << std::endl;
            std::cout << "-------------------------------------------------------------"  << std::endl;
        }

        //------------------------------------------------------------------
        // Particle Loop (Debug only for now)
        //------------------------------------------------------------------

        if (Debug){ // switch false with debug
            //for(int i = 0; i < bParticle -> GetEntriesFast(); ++i){
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
            
            //if (Debug) std::cout << " Good jet: " << goodjet << std::endl;

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
                hPr_e_eta -> Fill( lep-> Eta ,Event_Weight);
                hPr_e_Et -> Fill( TMath::Sqrt( TMath::Power(lep -> PT, 2) + TMath::Power(lep -> Mass, 2)), Event_Weight);
                
                aPr_e_eta -> FillWeighted(LepPass(lep), Event_Weight, lep->Eta);
                aPr_e_Et -> FillWeighted(LepPass(lep), Event_Weight, TMath::Sqrt( TMath::Power(lep -> PT, 2) + TMath::Power(lep -> Mass, 2)));

                //check if all 4e events seen 
                if(!LepPass(lep)) isalleseen = false;                
            }

            if( abs(lep->PID) == 12) {
                if (Debug) std::cout << " Nue E, Px, Py, Pz, P: " << Vec_Lepton.E() << " " << Vec_Lepton.Px() << " " << Vec_Lepton.Py() << " " << Vec_Lepton.Pz() << " " << Vec_Lepton.P()<<  std::endl;
                if (Debug) std::cout << " Nue theta, sintheta, costheta: "<<Vec_Lepton.Theta() <<  " " << TMath::Sin(Vec_Lepton.Theta()) << " " << TMath::Cos(Vec_Lepton.Theta()) <<std::endl;
                hPr_nue_eta -> Fill( lep-> Eta ,Event_Weight);
                hPr_nue_Et -> Fill( TMath::Sqrt( TMath::Power(lep -> PT, 2) + TMath::Power(lep -> Mass, 2)), Event_Weight);
            }

            //if it contains a muon, its definitely not a 4e event 
            if (abs(lep -> PID) == 13){
                is4e = false;
            }
        } // Lepton Loop

        
        //Lepton loop again but for 4e events only
        if (is4e){
            TLorentzVector Debug_MP;
            std::vector<GenParticle*> e_par_list;
            if (Debug) {
                std::cout << "This is a 4e event. Seen: " << isalleseen << std::endl;  
            }
            
            //into lepton for loop again, this time its only  4e events
            for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i){
                GenParticle* lep_e = (GenParticle*) bTruthLepton->At(i);
                TLorentzVector Vec_Lepton_e;
                Vec_Lepton_e.SetPtEtaPhiM(lep_e->PT,lep_e->Eta,lep_e->Phi,lep_e->Mass);

                if( abs(lep_e->PID) == 12) {
                    hEv_nue_eta_nocuts -> Fill( lep_e-> Eta ,Event_Weight);
                    hEv_nue_Et_nocuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);
                    if (isalleseen){
                        my_nu = lep_e;
                        //fill the electron neutrino comparison with missing energy
                        hEv_nue_eta_wicuts -> Fill( lep_e-> Eta ,Event_Weight);
                        hEv_nue_Et_wicuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);

                        hEv_nue_eta_pt_wicuts -> Fill(lep_e-> Eta,  lep_e-> PT);
                    }
                }
                
                if( abs(lep_e->PID) == 11){
                    hEv_e_eta_nocuts -> Fill( -(lep_e -> Eta) , Event_Weight);
                    hEv_e_Et_nocuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);

                    aEv_e_eta -> FillWeighted(LepPass(lep_e), Event_Weight, lep_e->Eta);
                    aEv_e_Et -> FillWeighted(LepPass(lep_e), Event_Weight, TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)));

                    if(LepPass(lep_e)){
                        hEv_e_eta_wicuts -> Fill( lep_e -> Eta , Event_Weight);
                        hEv_e_Et_wicuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);   
                        

                        e_par_list.push_back(lep_e);


                        Missing_Particle = Missing_Particle - Vec_Lepton_e;
                        FourLepton_Vector = FourLepton_Vector + Vec_Lepton_e;
                        if(Debug) {
                            //std::cout << "4elep " << i << " pT = " << Vec_Lepton_e.Pt() << " eta = " << Vec_Lepton_e.Eta() <<  std::endl;
                            //std::cout << " MPNow pT = " << Missing_Particle.Pt() << " eta = " << Missing_Particle.Eta() << std::endl;
                            //std::cout << "  4ee E, Px, Py, Pz: " << Vec_Lepton_e.E() << " " << Vec_Lepton_e.Px() << " " << Vec_Lepton_e.Py() << " " << Vec_Lepton_e.Pz() << std::endl;
                            //std::cout << "  MPe E, Px, Py, Pz: " << Missing_Particle.E() << " " << Missing_Particle.Px() << " " << Missing_Particle.Py() << " " << Missing_Particle.Pz() << std::endl;
                        }
    
                    }

                    hEv_e_eta_pt -> Fill(lep_e -> Eta, lep_e -> PT);
                }

            }

            hEx_EventCount -> Fill(1.5, Event_Weight);

            if (e_par_list.size() != 4) isalleseen = false; // so turns out there are some things that could possibly have 5 and i dont wanna deal with it
            
            if (isalleseen) {
                TLorentzVector my_nu_vec;
                
                my_nu_vec.SetPtEtaPhiM(my_nu->PT,my_nu->Eta,my_nu->Phi,my_nu->Mass);
                
                Debug_MP =   my_nu_vec - Missing_Particle;


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

                //logycut analysis // Andy's example Filling in Histograms

                for(int k = 0; k < nCuts; ++k){
                    if(TMath::Log10(H_Reco_Lst[2]) < cut_values.at(k)){   
                        h_varycut.at(k) -> Fill(FourLepton_Vector.M(), Event_Weight);
                    }
                }
               


                if (Debug){
                    std::cout << "MET pt: " << Missing_Particle.Pt()  << " eta: " << Missing_Particle.Eta() << std::endl;
                    //std::cout << "DebugMP pT = " << Debug_MP.Pt() << " eta = " << Debug_MP.Eta() << std::endl;
                    //std::cout << "DebugMP E, Px, Py, Pz: " << Debug_MP.E() << " " << Debug_MP.Px() << " " << Debug_MP.Py() << " " << Debug_MP.Pz() << std::endl;

                    std::cout << "E-Reco Q2: " << E_Reco_Lst[0] << " x: " << E_Reco_Lst[1] << " y: " << E_Reco_Lst[2] << std::endl;
                    std::cout << "H-Reco Q2: " << H_Reco_Lst[0] << " x: " << H_Reco_Lst[1] << " y: " << H_Reco_Lst[2] << std::endl;
                    std::cout << "Reco Higgs Mass?: " << FourLepton_Vector.M() << std::endl;
                    std::cout << "SumPz: " << sumEPz << std::endl;

                    std::cout << "ELIST CHECK " << e_par_list[0] -> Charge  << std::endl;

                    // TIME TO DO THE ZZ RECONSTRUCTION BIT HERE//
                    /*for(double d : ZZ_Reco(e_par_list).M()){
                        std::cout<<d<<std::endl;
                    }
                    */
                    
                    for(int oo = 0; oo < 2; ++oo){
                        //std::cout << e_par_list[oo] -> Charge << std::endl;
                        std::cout << Z_Reco_Lst[oo].M() << std::endl;
                    }
                    
                    
                    

                }

                hEx_EventCount -> Fill(2.5, Event_Weight);
                //if (true){ //filling stuff in histogram
                hEx_EventCount -> Fill(3.5, Event_Weight);

                
                hEv_MET_eta -> Fill(Missing_Particle.Eta() ,Event_Weight);
                hEv_MET_Et -> Fill(Missing_Particle.Pt() ,Event_Weight); //pt = et for neutrino
                
                hEv_MET_eta_pt -> Fill(Missing_Particle.Eta(), Missing_Particle.Pt());
                hEv_nue_MET_eta -> Fill(my_nu->Eta, Missing_Particle.Eta());
                hEv_nue_MET_Et -> Fill(my_nu -> PT, Missing_Particle.Pt());
                hEv_nue_MET_Phi -> Fill(my_nu -> Phi, Missing_Particle.Phi());

                hEv_debugMP_Pz_E -> Fill(Debug_MP.Pz(), Debug_MP.E());

                hEv_HReco_M -> Fill(FourLepton_Vector.M(), Event_Weight); //not the hadronic vector
                
                hEv_ZReco_M -> Fill(Z_Reco_Lst[0].M(), Event_Weight);
                hEv_ZstarReco_M -> Fill(Z_Reco_Lst[1].M(), Event_Weight);

                hEv_ZZReco_M -> Fill(Z_Reco_Lst[0].M(), Event_Weight);
                hEv_ZZReco_M -> Fill(Z_Reco_Lst[1].M(), Event_Weight);
    


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
                    
                    /*
                    } catch (...){
                        std::cout << "Caught exception while printing histogram" << std::endl;
                    }
                    */
                    
                
                //}
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
            aEv_H_eta ->  FillWeighted(isalleseen && is4e, Event_Weight, Higgs.Eta());
            
            //do the higgs selectoin events here 
            //as in acceptance of the higgs eta and pt
	    }
    } // Loop over all events
}

