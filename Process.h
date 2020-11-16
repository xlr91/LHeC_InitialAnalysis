#include "TFile.h"
#include "TH2.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"
#include "TStyle.h"
#include "TVector.h"
#include "TError.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TVector3.h"
#include "TEfficiency.h"
#include "TLine.h"

#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"
#include "TRandom3.h"
#include "TClonesArray.h"

#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

// Plots

TClonesArray * bEvent;
TClonesArray * bTruthLepton;
TClonesArray * bTruthWZ;
TClonesArray * bJet;
TClonesArray * bParticle;




// Output
TFile * OutputFile;

TH1D * hEx_EventCount;
TH1D * hEx_WeightCount;

TH1D * hEx_Z_Pt;
TH1D * hEx_Lepton_Pt;
TH1D * hEx_Jet_Pt;
TH1D * hEx_ZZ_Mass;

//Particle level graphs
TH1D * hPr_e_eta;
TH1D * hPr_nue_eta;
TH1D * hPr_Jet_eta;
TH1D * hPr_Z_eta;
TH1D * hPr_H_eta;

TH1D * hPr_e_Et;
TH1D * hPr_nue_Et;
TH1D * hPr_Jet_Et;
TH1D * hPr_Z_Et;

TEfficiency * aPr_e_eta;
TEfficiency * aPr_e_Et;


//4e event level graphs
TH1D * hEv_e_eta_nocuts;
TH1D * hEv_e_eta_wicuts;
TH1D * hEv_e_Et_nocuts;
TH1D * hEv_e_Et_wicuts;
TH1D * hEv_nue_eta_nocuts;
TH1D * hEv_nue_Et_nocuts;
TH1D * hEv_nue_eta_wicuts;
TH1D * hEv_nue_Et_wicuts;
TH1D * hEv_MET_eta;
TH1D * hEv_MET_Et;

TH2D * hEv_e_eta_pt;

TH2D * hEv_nue_eta_pt_wicuts;
TH2D * hEv_MET_eta_pt;

TEfficiency * aEv_H_eta;

TEfficiency * aEv_e_eta;
TEfficiency * aEv_e_Et;
            


ExRootTreeReader * InitReader(const TString FilePath);

void Process(ExRootTreeReader * treeReader);

void ClearBranches();

int main(int argc, char* argv[]);

bool LepPass(GenParticle* lep_b);
//double jetlepdR (TLorentzVector jet_v, TLorentzVector lep_v)
