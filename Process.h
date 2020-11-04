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


// Output
TFile * OutputFile;

TH1D * h_EventCount;
TH1D * h_WeightCount;

TH1D * h_Z_Pt;
TH1D * h_Lepton_Pt;
TH1D * h_Jet_Pt;
TH1D * h_ZZ_Mass;

//my graphs
TH1D * h_e_eta;
TH1D * h_nue_eta;
TH1D * h_Jet_eta;
TH1D * h_Z_eta;
TH1D * h_H_eta;

TH1D * h_e_Et;
TH1D * h_nue_Et;
TH1D * h_Jet_Et;
TH1D * h_Z_Et;

//TEfficiency for acceptance
TEfficiency * t_e_eta;
TEfficiency * t_e_Et;


ExRootTreeReader * InitReader(const TString FilePath);

void Process(ExRootTreeReader * treeReader);

void ClearBranches();

int main(int argc, char* argv[]);

bool LepPass(GenParticle* lep_b);

