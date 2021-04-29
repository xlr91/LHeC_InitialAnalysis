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

TH2D * hEv_nue_MET_eta;
TH2D * hEv_nue_MET_Et;
TH2D * hEv_nue_MET_Phi;

TH2D * hEv_debugMP_Pz_E;

TEfficiency * aEv_H_eta;

TEfficiency * aEv_e_eta;
TEfficiency * aEv_e_Et;

TH1D * hEv_HReco_M;
TH1D * hEv_ZReco_M;
TH1D * hEv_ZstarReco_M;
TH1D * hEv_ZZReco_M;
TH1D * hEv_ZCurr_M;

TH1D * hEv_jet_eta;
TH1D * hEv_H_eta;



TH2D * hEvR_recoQ2_elec_hadr;
TH2D * hEvR_recox_elec_hadr;
TH2D * hEvR_recoy_elec_hadr;

TF1 * hEvR_recoQ2_elec_hadr_fit;



TH1D * hEvR_ereco_Q2;
TH1D * hEvR_ereco_x;
TH1D * hEvR_ereco_y;
TH1D * hEvR_ereco_y_true;

TH1D * hEvR_hreco_Q2;
TH1D * hEvR_hreco_x;
TH1D * hEvR_hreco_y;
TH1D * hEvR_hreco_y_true;
TH1D * hEvR_EPz;

TH2D * hEvR_hreco_x_Q2;


TH1D * hEvC_Zstar;
TH1D * hEvC_Logy;

TH1D * hEvS_e_pt;
TH1D * hEvS_e_pt_S;

TH1D * hEvS_e_E;
TH1D * hEvS_e_E_S;

TH1D * hEv_HReco_M_S;
TH1D * hEv_ZReco_M_S;
TH1D * hEv_ZstarReco_M_S;

std::vector<TH1D*> h_varycut;
std::vector<double> cut_values;
int nCuts;
double min_cut; 


            


ExRootTreeReader * InitReader(const TString FilePath);

void Process(ExRootTreeReader * treeReader, TString Ident);

void ClearBranches();

int main(int argc, char* argv[]);

//bool LepPass(GenParticle* lep_b);
bool LepPass(GenParticle* lep_b, double ptsmear);

std::vector<double> Electron_Reco(TLorentzVector scat);
std::vector<double> Hadron_Reco(TLorentzVector had);

//std::vector<TLorentzVector> ZZ_Reco(std::vector<GenParticle*> test);
std::tuple<std::vector<TLorentzVector>,int, std::vector<TLorentzVector>,int> ZZ_Reco(std::vector<GenParticle*> test, std::vector<TLorentzVector> e_vec_smear);

double ptSmear(TRandom* gR, GenParticle* lep_s);
double ESmear(TRandom* gR, GenParticle* lep_s);
//double jetlepdR (TLorentzVector jet_v, TLorentzVector lep_v)
