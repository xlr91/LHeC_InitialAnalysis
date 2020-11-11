
#include "Process.h"
#include "TLine.h"

bool Debug = true;

bool LepPass(GenParticle* lep_b){
    //Parameters
    double etamin = -5 ;
    double etamax = 5;
    double ptmin = 5;

    //selection criteria
    if( etamin < lep_b -> Eta && lep_b -> Eta<  etamax && ptmin < lep_b -> PT){
        return true;
    }
    return false;
}  

int main(int argc, char* argv[]) {

    // Input Delphes File

    const TString InputFile = argv[1];
    const TString OutputFileName = argv[2];

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Running Process with my edit"  << std::endl;
    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "InputFile = " << InputFile << std::endl;
    std::cout << "OutputFileName = " << OutputFileName << std::endl;
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

    hEx_EventCount = new TH1D("hEx_EventCount","Event Classifications ; Event type; Number of Events",4,0,4);
    hEx_WeightCount = new TH1D("hEx_WeightCount","",1,0,1);

    hEx_Lepton_Pt = new TH1D("hEx_Lepton_Pt","Charged Lepton Events vs pT; Charged Lepton p_{T} [GeV]; Number of Particles / 5 GeV",200,0.0,1000.0);
    hEx_Z_Pt = new TH1D("hEx_Z_Pt","Z Boson Events vs pT; Z p_{T} [GeV]; Number of Particles / 5 GeV",200,0.0,1000.0);
    hEx_Jet_Pt = new TH1D("hEx_Jet_Pt","Jet Events vs pT; Jet p_{T} [GeV]; Number of Particles / 5 GeV",200,0.0,1000.0);
    hEx_ZZ_Mass = new TH1D("hEx_ZZ_Mass","Events vs mass of ZZ*; ZZ Mass [GeV]; Number of Particles / 2 GeV",125,0.0,250.0);


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


    hEv_e_eta_nocuts = new TH1D("hEv_e_eta_nocuts","Electron (4e, no cuts) particles vs pseudorapidity ; Electron #eta; Number of Particles", 50, -7.0, 3.0);
    hEv_e_eta_wicuts = new TH1D("hEv_e_eta_wicuts","Electron (4e, with cuts) particles vs pseudorapidity; Electron #eta; Number of Particles", 50, -7.0, 3.0);
    hEv_e_Et_nocuts = new TH1D("hEv_e_Et_nocuts","Electron (4e, no cuts) particles vs transverse energy; Electron E_{T}; Number of Particles", 50, 0, 150);
    hEv_e_Et_wicuts = new TH1D("hEv_e_Et_wicuts","Electron (4e, with cuts) particles vs transverse energy; Electron E_{T} ; Number of Particles", 50, 0, 150);
    hEv_nue_eta = new TH1D("hEv_nue_eta","Electron Neutrino (4e, no cuts) particles vs pseudorapidity; Electron neutrino #eta; Number of Particles", 50, -5.0, 5.0);
    hEv_nue_Et = new TH1D("hEv_nue_Et","Electron Neutrino (4e, no cuts) particles vs transverse energy; Electron Neutrino E_{T}; Number of Particles", 50, 0, 250);
    hEv_MET_eta = new TH1D("hEv_MET_eta","Missing Energy Particle (with all cuts) vs pseudorapidity; MET #eta; Number of Particles", 50, -5.0, 5.0);
    hEv_MET_Et = new TH1D("hEv_MET_Et","Missing Energy Particle (with all cuts) vs transverse energy; MET E_{T}; Number of Particles", 50, 0, 250);
    aEv_e_eta = new TEfficiency("aEv_e_eta","Acceptance of Electron vs Eta (4e events); Eta; Acceptance",100,-10,10);
    aEv_e_Et = new TEfficiency("aEv_e_Et","Acceptance of Electron vs Transverse Energy (4e events); Et; Acceptance", 50 , 0, 150);

    hEv_e_eta_pt = new TH2D("hEv_e_eta_pt","", 50, -7.0, 3.0, 50, 0., 150.);

    //------------------------------------

    // Run the selection
    Process(reader);


    //Writes selection to the output file
    ///std::cout << "Events in EventCount: " << hEx_EventCount->GetEntries() << std::endl;
    std::cout << "Events in EventCount: " << hEx_EventCount->GetBinContent(1) << std::endl;
    std::cout << "Events that are 4e  : " << hEx_EventCount->GetBinContent(2) << std::endl;
    std::cout << "Events seen thats 4e: " << hEx_EventCount->GetBinContent(3) << std::endl;
    std::cout << "Events seen 4e GJet : " << hEx_EventCount->GetBinContent(4) << std::endl;

    

    std::cout << "Write to file..." << std::endl;


    //draw lines here
    /*
        TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,1680,1050);
        h_Jet_eta -> Draw();
        TLine *line1 = new TLine(-5,0,-5,2000);
        line1->SetLineColor(kRed);
        line1->Draw("same");
        c1 -> Print("testing.png");
    */

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
    hEv_nue_eta -> Write();
    hEv_nue_Et -> Write();
    hEv_MET_eta -> Write();
    hEv_MET_Et -> Write();

    aEv_e_eta -> Write(); 
    aEv_e_Et -> Write(); 

    hEv_e_eta_pt -> Write();


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

void Process(ExRootTreeReader * treeReader) {

    // Get pointers to branches used in this analysis
    bEvent = treeReader->UseBranch("Event");
    bJet = treeReader->UseBranch("GenJet");
    bTruthLepton = treeReader->UseBranch("TruthLeptonParticles");
    bTruthWZ = treeReader->UseBranch("TruthWZHParticles");


    bParticle = treeReader -> UseBranch("Particle");


    Long64_t numberOfEntries = treeReader->GetEntries();

    bool is4e;
    bool isall4eseen;
    bool goodjetevent;
    bool goodjet;
    int goodjetnum;

    if (Debug) numberOfEntries = 1000;

    int nSelected = 0;

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Input: " << numberOfEntries << " events to process" << std::endl;

    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        HepMCEvent * event = (HepMCEvent*) bEvent->At(0); // what does this mean
    	const float Event_Weight = event->Weight; //what is weight

        hEx_EventCount->Fill(0.5);
        hEx_WeightCount->Fill(0.5,Event_Weight);

        if( (entry > 0 && entry%1000 == 0) || Debug) {
            std::cout << "-------------------------------------------------------------"  << std::endl;
            std::cout << "Processing Event Number =  " << entry  << std::endl;
            std::cout << "-------------------------------------------------------------"  << std::endl;
        }

        //------------------------------------------------------------------
        // Particle Loop
        //------------------------------------------------------------------

        if (Debug){
            for(int i = 0; i < bParticle -> GetEntriesFast(); ++i){
                GenParticle* p_Particle = (GenParticle*) bParticle->At(i);
                std::cout << "Particle " << i << " pT = " << p_Particle->PT << " eta = " << p_Particle->Eta << " phi = " << p_Particle->Phi << " PID = " << p_Particle-> PID << " mass = " << p_Particle->Mass
                    << " Mother: " << p_Particle->M1 << " Daughter: " << p_Particle->D1 << " Status: " << p_Particle -> Status <<  std::endl;
            }
        }

        //------------------------------------------------------------------
        // Jet Loop
        //------------------------------------------------------------------
        goodjetnum = 0;
        goodjetevent = true;
        TLorentzVector Vec_Jet_Good;

        for(int i = 0; i < bJet->GetEntriesFast(); ++i) {

            Jet* jet = (Jet*) bJet->At(i);
            
            TLorentzVector Vec_Jet;

            Vec_Jet.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);

		    hEx_Jet_Pt->Fill( Vec_Jet.Pt(), Event_Weight ); // event weight is required apparently
            hPr_Jet_eta -> Fill( jet-> Eta ,Event_Weight);
            hPr_Jet_Et -> Fill( TMath::Sqrt( TMath::Power(jet -> PT, 2) + TMath::Power(jet -> Mass, 2)), Event_Weight);

            //do the cuts here
            goodjet = true;
            if (jet->Mass < 0.001) goodjet = false;
            if (goodjet){ 
                goodjetnum = goodjetnum + goodjet;
                Vec_Jet_Good = Vec_Jet;
                }

            if(Debug) std::cout << "Jet " << i << " pT = " << jet->PT << " eta = " << jet->Eta << " phi = " << jet->Phi << " mass = " << jet->Mass << " flavour = " << jet->Flavor  
                << " GoodJet: " << goodjet << std::endl;
        } // Jet Loop
        
        if (goodjetnum != 1) goodjetevent = false;

        //------------------------------------------------------------------
        // Lepton Loop
        //------------------------------------------------------------------

        is4e = true; //flag used to check if the lepton event satisfies the 4e subdecay channel
        isall4eseen = true;
        TLorentzVector Vec_Lepton_f;


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
                if( !LepPass(lep)) isall4eseen = false;
                Vec_Lepton_f = Vec_Lepton_f + Vec_Lepton; //note: does not take into account possiblity of 6 electrons or more yet
            }

            if( abs(lep->PID) == 12) {
                hPr_nue_eta -> Fill( lep-> Eta ,Event_Weight);
                hPr_nue_Et -> Fill( TMath::Sqrt( TMath::Power(lep -> PT, 2) + TMath::Power(lep -> Mass, 2)), Event_Weight);
            }

            //4e event selection here
            if (abs(lep -> PID) == 13){
                is4e = false;
            }
        } // Lepton Loop

        Vec_Lepton_f = Vec_Lepton_f + Vec_Jet_Good;
        if (Debug) std::cout << "fourvector addition test " << Vec_Lepton_f.Pt() << std::endl;
        
        //Lepton loop again but for 4e events only
        if (is4e){
            if (Debug) {
                std::cout << "This is a 4e event. Seen: " << isall4eseen << std::endl;
                std::cout << "Number of good jets: " << goodjetnum << std::endl;
            }
            
            for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i){
                GenParticle* lep_e = (GenParticle*) bTruthLepton->At(i);

                if( abs(lep_e->PID) == 12) {
                    hEv_nue_eta -> Fill( lep_e-> Eta ,Event_Weight);
                    hEv_nue_Et -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);
                }

                if( abs(lep_e->PID) == 11){
                    hEv_e_eta_nocuts -> Fill( lep_e -> Eta , Event_Weight);
                    hEv_e_Et_nocuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);

                    aEv_e_eta -> FillWeighted(LepPass(lep_e), Event_Weight, lep_e->Eta);
                    aEv_e_Et -> FillWeighted(LepPass(lep_e), Event_Weight, TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)));

                    if(LepPass(lep_e)){
                        hEv_e_eta_wicuts -> Fill( lep_e -> Eta , Event_Weight);
                        hEv_e_Et_wicuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);
                        
                    }

                    hEv_e_eta_pt -> Fill(lep_e -> Eta, lep_e -> PT);
                }

            }

            hEx_EventCount -> Fill(1.5);
            if (isall4eseen) {
                hEx_EventCount -> Fill(2.5);
                if (goodjetevent){ 
                    hEx_EventCount -> Fill(3.5);}
                    hEv_MET_eta -> Fill(Vec_Lepton_f.Eta() ,Event_Weight);
                    hEv_MET_Et -> Fill(Vec_Lepton_f.Pt() ,Event_Weight); //pt = et for neutrino
                
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
	    }


    } // Loop over all events


}

