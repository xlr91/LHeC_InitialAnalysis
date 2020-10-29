
#include "Process.h"

bool Debug = true; 

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

    h_EventCount = new TH1D("h_EventCount","",1,0,1);
    h_WeightCount = new TH1D("h_WeightCount","",1,0,1);

    h_Lepton_Pt = new TH1D("h_Lepton_Pt","; Charged Lepton p_{T} [GeV]; Events / 5 GeV",200,0.0,1000.0);
    h_Z_Pt = new TH1D("h_Z_Pt","; Z p_{T} [GeV]; Events / 5 GeV",200,0.0,1000.0);
    h_Jet_Pt = new TH1D("h_Jet_Pt","; Jet p_{T} [GeV]; Events / 5 GeV",200,0.0,1000.0);
    h_ZZ_Mass = new TH1D("h_ZZ_Mass","; ZZ Mass [GeV]; Events / 2 GeV",125,0.0,250.0);


    //------------------------------------

    // Run the selection
    Process(reader);
    //Writes selection to the output file
    std::cout << "Events in EventCount: " << h_EventCount->GetEntries() << std::endl;

    std::cout << "Write to file..." << std::endl;

    OutputFile->cd();

    h_EventCount->Write();
    h_WeightCount->Write();

    h_Z_Pt->Write();
    h_Lepton_Pt->Write();
    h_Jet_Pt->Write();
    h_ZZ_Mass->Write();

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

    int nSelected = 0;

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Input: " << numberOfEntries << " events to process" << std::endl;

    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        HepMCEvent * event = (HepMCEvent*) bEvent->At(0);
    	const float Event_Weight = event->Weight;

        h_EventCount->Fill(0.5);
        h_WeightCount->Fill(0.5,Event_Weight);

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

		h_Jet_Pt->Fill( Vec_Jet.Pt(), Event_Weight );

        } // Jet Loop


        //------------------------------------------------------------------
        // Lepton Loop
        //------------------------------------------------------------------
        for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) {

            GenParticle* lep = (GenParticle*) bTruthLepton->At(i);

            if(Debug) std::cout << "Lepton " << i << " pT = " << lep->PT << " eta = " << lep->Eta << " phi = " << lep->Phi << " mass = " << lep->Mass << 
std::endl;

            
            TLorentzVector Vec_Lepton;
            Vec_Lepton.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);
            

		// Look for electrons or muons
		if( abs(lep->PID) == 11 || abs(lep->PID) == 13 ) {

	            h_Lepton_Pt->Fill( Vec_Lepton.Pt(), Event_Weight );  

		}

        } // Lepton Loop

	std::vector<TLorentzVector> list_Zboson;

        //------------------------------------------------------------------
        // W/Z/H Boson Loop
        //------------------------------------------------------------------
        for(int i = 0; i < bTruthWZ->GetEntriesFast(); ++i) {
            
            GenParticle* v = (GenParticle*) bTruthWZ->At(i);
            
            if(Debug) std::cout << "Lepton " << i << " pT = " << v->PT << " eta = " << v->Eta << " phi = " << v->Phi << " mass = " << v->Mass <<
std::endl;
        
        
            TLorentzVector Vec_V;
            Vec_V.SetPtEtaPhiM(v->PT,v->Eta,v->Phi,v->Mass);
        
    		// Look for Z bosons (ID = 23)
            	if( v->PID == 23 ) {
    
                    h_Z_Pt->Fill( Vec_V.Pt(), Event_Weight );
    		
			list_Zboson.push_back(Vec_V);
             }
            
        } // W/Z/H Loop

	if( list_Zboson.size() > 1 ) {

		TLorentzVector Higgs = list_Zboson.at(0) + list_Zboson.at(1);

              h_ZZ_Mass->Fill( Higgs.M() , Event_Weight );
	}


    } // Loop over all events


}
