#include "Process.h"

//fix the labels of the plots

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

    //Example histograms
    hEx_EventCount = new TH1D("hEx_EventCount","Event Classifications ; Event type; Number of Events",4,0,4);
    hEx_WeightCount = new TH1D("hEx_WeightCount","",1,0,1);

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
    hEv_e_eta_nocuts = new TH1D("hEv_e_eta_nocuts","Electron (4e, no cuts) particles vs pseudorapidity ; Electron #eta; Number of Particles", 50, -7.0, 3.0);
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

    hEv_e_eta_pt = new TH2D("hEv_e_eta_pt","", 50, -7.0, 3.0, 50, 0., 150.);
    hEv_nue_eta_pt_wicuts = new TH2D("hEv_nue_eta_pt_wicuts","", 50, -10.0, 10.0, 50, 0., 250.);
    hEv_MET_eta_pt = new TH2D("hEv_MET_eta_pt","", 50, -10.0, 10.0, 50, 0., 250.);
    hEv_nue_MET_Phi = new TH2D("hEv_nue_MET_Phi","", 50, -4.0, 4.0, 50, -4., 4.);

    hEv_nue_MET_eta = new TH2D("hEv_nue_MET_eta","", 100, -5.0, 10.0, 100, -5.0, 10.0);
    hEv_nue_MET_Et = new TH2D("hEv_nue_MET_Et","", 50, 0, 250, 50, 0, 250);


    hEv_debugMP_Pz_E = new TH2D("hEv_debugMP_Pz_E","", 100, -2000, 0, 100, 0, 2000);


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
    if (Debug) numberOfEntries = 1000;

    bool is4e; //is the event a 4e event?
    bool isalleseen; //can all electrons in the event be seen?
    bool goodjet; //is this jet a good jet? (Does its mass *not* correspond to a lepton)

    int nSelected = 0;

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Input: " << numberOfEntries << " events to process" << std::endl;

    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        HepMCEvent * event = (HepMCEvent*) bEvent->At(0); 
    	const float Event_Weight = event->Weight;

        hEx_EventCount->Fill(0.5);
        hEx_WeightCount->Fill(0.5,Event_Weight);

        is4e = true;
        isalleseen = true;
        TLorentzVector Vec_Lepton_f;
        TLorentzVector Missing_Particle;
        Missing_Particle.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        GenParticle* my_nu;

        if( (entry > 0 && entry%1000 == 0) || Debug) {
            std::cout << "-------------------------------------------------------------"  << std::endl;
            std::cout << "Processing Event Number =  " << entry  << std::endl;
            std::cout << "-------------------------------------------------------------"  << std::endl;
        }

        //------------------------------------------------------------------
        // Particle Loop (Debug only for now)
        //------------------------------------------------------------------

        if (false){ // switch false with debug
            for(int i = 0; i < bParticle -> GetEntriesFast(); ++i){
                GenParticle* p_Particle = (GenParticle*) bParticle->At(i);
                std::cout << "Particle " << i << " pT = " << p_Particle->PT << " eta = " << p_Particle->Eta << " phi = " << p_Particle->Phi << " PID = " << p_Particle-> PID << " mass = " << p_Particle->Mass
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
            
            if (Debug) std::cout << " Good jet: " << goodjet << std::endl;

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
                if (Debug) std::cout << " Nue E, Px, Py, Pz: " << Vec_Lepton.E() << " " << Vec_Lepton.Px() << " " << Vec_Lepton.Py() << " " << Vec_Lepton.Pz() << std::endl;
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
                    hEv_e_eta_nocuts -> Fill( lep_e -> Eta , Event_Weight);
                    hEv_e_Et_nocuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);

                    aEv_e_eta -> FillWeighted(LepPass(lep_e), Event_Weight, lep_e->Eta);
                    aEv_e_Et -> FillWeighted(LepPass(lep_e), Event_Weight, TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)));

                    if(LepPass(lep_e)){
                        hEv_e_eta_wicuts -> Fill( lep_e -> Eta , Event_Weight);
                        hEv_e_Et_wicuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);   
                        
                        Missing_Particle = Missing_Particle - Vec_Lepton_e;
                        if(Debug) {
                            std::cout << "4elep " << i << " pT = " << Vec_Lepton_e.Pt() << " eta = " << Vec_Lepton_e.Eta() <<  std::endl;

                            std::cout << " MPNow pT = " << Missing_Particle.Pt() << " eta = " << Missing_Particle.Eta() << std::endl;
                            std::cout << "  4ee E, Px, Py, Pz: " << Vec_Lepton_e.E() << " " << Vec_Lepton_e.Px() << " " << Vec_Lepton_e.Py() << " " << Vec_Lepton_e.Pz() << std::endl;
                            std::cout << "  MPe E, Px, Py, Pz: " << Missing_Particle.E() << " " << Missing_Particle.Px() << " " << Missing_Particle.Py() << " " << Missing_Particle.Pz() << std::endl;
                        }
    
                    }

                    hEv_e_eta_pt -> Fill(lep_e -> Eta, lep_e -> PT);
                }

            }

            TLorentzVector my_nu_vec;
            my_nu_vec.SetPtEtaPhiM(my_nu->PT,my_nu->Eta,my_nu->Phi,my_nu->Mass);
            Debug_MP =   my_nu_vec - Missing_Particle;

            if (Debug){
                std::cout << "MET pt: " << Missing_Particle.Pt()  << " eta: " << Missing_Particle.Eta() << std::endl;
                


                std::cout << "DebugMP pT = " << Debug_MP.Pt() << " eta = " << Debug_MP.Eta() << std::endl;
                std::cout << "DebugMP E, Px, Py, Pz: " << Debug_MP.E() << " " << Debug_MP.Px() << " " << Debug_MP.Py() << " " << Debug_MP.Pz() << std::endl;
            }

            hEx_EventCount -> Fill(1.5);
            
            if (isalleseen) {
                hEx_EventCount -> Fill(2.5);

                if (true){ 
                    hEx_EventCount -> Fill(3.5);

                    
                    hEv_MET_eta -> Fill(Missing_Particle.Eta() ,Event_Weight);
                    hEv_MET_Et -> Fill(Missing_Particle.Pt() ,Event_Weight); //pt = et for neutrino
                    
                    hEv_MET_eta_pt -> Fill(Missing_Particle.Eta(), Missing_Particle.Pt());
                    hEv_nue_MET_eta -> Fill(my_nu->Eta, Missing_Particle.Eta());
                    hEv_nue_MET_Et -> Fill(my_nu -> PT, Missing_Particle.Pt());
                    hEv_nue_MET_Phi -> Fill(my_nu -> Phi, Missing_Particle.Phi());

                    
                    hEv_debugMP_Pz_E -> Fill(Debug_MP.Pz(), Debug_MP.E());
                    
                
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
            aEv_H_eta ->  FillWeighted(isalleseen && is4e, Event_Weight, Higgs.Eta());
            

            //do the higgs selectoin events here 
            //as in acceptance of the higgs eta and pt
	    }


    } // Loop over all events


}

