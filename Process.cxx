
#include "Process.h"
#include "TLine.h"

bool Debug = false; 

bool LepPass(GenParticle* lep_b){
    //Parameters
    double etamin = -5 ;
    double etamax = 5;
    double ptmin = 15;

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

    hEx_EventCount = new TH1D("hEx_EventCount","",3,0,3);
    hEx_WeightCount = new TH1D("hEx_WeightCount","",1,0,1);

    hEx_Lepton_Pt = new TH1D("hEx_Lepton_Pt","; Charged Lepton p_{T} [GeV]; Events / 5 GeV",200,0.0,1000.0);
    hEx_Z_Pt = new TH1D("hEx_Z_Pt","; Z p_{T} [GeV]; Events / 5 GeV",200,0.0,1000.0);
    hEx_Jet_Pt = new TH1D("hEx_Jet_Pt","; Jet p_{T} [GeV]; Events / 5 GeV",200,0.0,1000.0);
    hEx_ZZ_Mass = new TH1D("hEx_ZZ_Mass","; ZZ Mass [GeV]; Events / 2 GeV",125,0.0,250.0);


    hPr_Jet_eta = new TH1D("hPr_Jet_eta","; Jet eta ; Events", 50, -10.0, 0.0);
    hPr_e_eta = new TH1D("hPr_e_eta","; Electron eta ; Events", 50, -7.0, 3.0);
    hPr_nue_eta = new TH1D("hPr_nue_eta","; Electron neutrino eta ; Events", 50, -5.0, 5.0);
    hPr_Z_eta = new TH1D("hPr_Z_eta","; Boson eta ; Events", 50, -10.0, 0.0);
    hPr_H_eta = new TH1D("hPr_H_eta","; Higgs eta ; Events", 100, -10.0, 10.0);

    hPr_Jet_Et = new TH1D("hPr_Jet_Et","; Jet Transverse Energy ; Events", 50, 0, 200);
    hPr_e_Et = new TH1D("hPr_e_Et","; Electron Transverse Energy; Events", 50, 0, 150);
    hPr_nue_Et = new TH1D("hPr_nue_Et","; Electron Neutrino Transverse Energy; Events", 50, 0, 250);
    hPr_Z_Et = new TH1D("hPr_Z_Et","; Boson Transverse Energy ; Events", 50, 0, 250);

    tPr_e_eta = new TEfficiency("tPr_e_eta","Acceptance of Electron vs Eta ; Eta; Acceptance",100,-10,10);
    tPr_e_Et = new TEfficiency("tPr_e_Et","Acceptance of Electron vs Transverse Energy ; Et; Acceptance", 50 , 0, 150);


    hEv_e_eta_nocuts = new TH1D("hEv_e_eta_nocuts","; Electron eta (4e events, no cuts); Events", 50, -7.0, 3.0);
    hEv_e_eta_wicuts = new TH1D("hEv_e_eta_wicuts","; Electron eta (4e events, with cuts); Events", 50, -7.0, 3.0);
    hEv_nue_eta = new TH1D("hEv_nue_eta","; Electron neutrino eta (4e events, no cuts); Events", 50, -5.0, 5.0);
    hEv_e_Et_nocuts = new TH1D("hEv_e_Et_nocuts","; Electron Transverse Energy (4e events, no cuts); Events", 50, 0, 150);
    hEv_e_Et_wicuts = new TH1D("hEv_e_Et_wicuts","; Electron Transverse Energy (4e events, with cuts); Events", 50, 0, 150);
    hEv_nue_Et = new TH1D("hEv_nue_Et","; Electron Neutrino Transverse Energy (4e events, no cuts); Events", 50, 0, 250);
    tEv_e_eta = new TEfficiency("tEv_e_eta","Acceptance of Electron vs Eta (4e events); Eta; Acceptance",100,-10,10);
    tEv_e_Et = new TEfficiency("tEv_e_Et","Acceptance of Electron vs Transverse Energy (4e events); Et; Acceptance", 50 , 0, 150);



    //------------------------------------

    // Run the selection
    Process(reader);


    //Writes selection to the output file
    ///std::cout << "Events in EventCount: " << hEx_EventCount->GetEntries() << std::endl;
    std::cout << "Events in EventCount: " << hEx_EventCount->GetBinContent(1) << std::endl;
    std::cout << "Events that are 4e  : " << hEx_EventCount->GetBinContent(2) << std::endl;
    std::cout << "Events seen thats 4e: " << hEx_EventCount->GetBinContent(3) << std::endl;

    

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

    tPr_e_eta -> Write(); 
    tPr_e_Et -> Write(); 


    OutputFile->cd("4eEventLevel");

    hEv_e_eta_nocuts -> Write();
    hEv_e_eta_wicuts -> Write();
    hEv_nue_eta -> Write();
    hEv_e_Et_nocuts -> Write();
    hEv_e_Et_wicuts -> Write();
    hEv_nue_Et -> Write();

    tEv_e_eta -> Write(); 
    tEv_e_Et -> Write(); 


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

    Long64_t numberOfEntries = treeReader->GetEntries();

    bool is4e;
    bool isall4eseen;

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
        // Jet Loop
        //------------------------------------------------------------------

        for(int i = 0; i < bJet->GetEntriesFast(); ++i) {

            Jet* jet = (Jet*) bJet->At(i);

            if(Debug) std::cout << "Jet " << i << " pT = " << jet->PT << " eta = " << jet->Eta << " phi = " << jet->Phi << " mass = " << jet->Mass << " flavour = " << jet->Flavor << std::endl;

            TLorentzVector Vec_Jet;
            Vec_Jet.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);

		    hEx_Jet_Pt->Fill( Vec_Jet.Pt(), Event_Weight ); // event weight is required apparently
            hPr_Jet_eta -> Fill( jet-> Eta ,Event_Weight);
            hPr_Jet_Et -> Fill( TMath::Sqrt( TMath::Power(jet -> PT, 2) + TMath::Power(jet -> Mass, 2)), Event_Weight);

        } // Jet Loop


        //------------------------------------------------------------------
        // Lepton Loop
        //------------------------------------------------------------------

        is4e = true; //flag used to check if the lepton event satisfies the 4e subdecay channel
        isall4eseen = true;
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
                
                tPr_e_eta -> FillWeighted(LepPass(lep), Event_Weight, lep->Eta);
                tPr_e_Et -> FillWeighted(LepPass(lep), Event_Weight, TMath::Sqrt( TMath::Power(lep -> PT, 2) + TMath::Power(lep -> Mass, 2)));
                ///if (Debug) std::cout << "Selection result: " << LepPass(lep) << std::endl;

                //check if all 4e events seen 
                if( !LepPass(lep)) isall4eseen = false;

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

        //Lepton loop again but for 4e events only
        if (is4e){
            if (Debug) {
                std::cout << "This is a 4e event. Seen: " << isall4eseen << std::endl;
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

                    tEv_e_eta -> FillWeighted(LepPass(lep_e), Event_Weight, lep_e->Eta);
                    tEv_e_Et -> FillWeighted(LepPass(lep_e), Event_Weight, TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)));

                    if(LepPass(lep_e)){
                        hEv_e_eta_wicuts -> Fill( lep_e -> Eta , Event_Weight);
                        hEv_e_Et_wicuts -> Fill( TMath::Sqrt( TMath::Power(lep_e -> PT, 2) + TMath::Power(lep_e -> Mass, 2)), Event_Weight);
                    }
                }
            }

            hEx_EventCount -> Fill(1.5);
            if (isall4eseen) hEx_EventCount -> Fill(2.5);
            
            /*
            for electrons 
                fill eta
                fill et
                //do we also want one with the cuts?
                fill acceptances

                if isall4eseen
                    more plots here if you want
            */

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
