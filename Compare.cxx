#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <TCanvas.h>
#include "TF1.h"
#include "TEfficiency.h"

void Compare(){
    TString filename1 = "signal.root";
    TString filename2 = "bakgnd.root";

    TFile *f1 = TFile::Open(filename1);
    TFile *f2 = TFile::Open(filename2);
    TFile *OutputFile = new TFile("Compare.root", "recreate");
    OutputFile -> cd();
    //TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,1680,1050); //original
    TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,640,480); //for quicklook
    gStyle -> SetOptStat(0);

    TH1D *h1, *h2;
    TLegend *legend;
    Double_t histmax;

    
    std::vector<TString> h1_names;

    h1_names.push_back("4eEventLevel/hEv_HReco_M");
    h1_names.push_back("4eEventLevel/hEv_ZReco_M");
    h1_names.push_back("4eEventLevel/hEv_ZstarReco_M");
    h1_names.push_back("4eEventLevel/hEv_MET_Et");
    h1_names.push_back("4eEventLevel/hEv_MET_eta");
    h1_names.push_back("4eEventLevel/hEv_e_eta_nocuts");
    h1_names.push_back("4eEventLevel/hEv_e_Et_nocuts");

    h1_names.push_back("ParticleLevel/hPr_Jet_eta");
    h1_names.push_back("ParticleLevel/hPr_Jet_Et");
    

    h1_names.push_back("4eEventLevel/KinematicReco/hEvR_hreco_x");
    h1_names.push_back("4eEventLevel/KinematicReco/hEvR_hreco_y");
    h1_names.push_back("4eEventLevel/KinematicReco/hEvR_hreco_Q2");


    //TH1 loop

    for(Int_t k = 0; k < h1_names.size() ; k++){
        h1 = (TH1D*) f1 -> Get(h1_names[k]);
        h2 = (TH1D*) f2 -> Get(h1_names[k]);    
        
        histmax = h1 -> GetMaximum();
        if(histmax < h2 -> GetMaximum()) histmax = h2 -> GetMaximum();

        h1 -> SetAxisRange(0, histmax*1.1, "Y");
        h1 -> SetLineColor(2);

        h1 -> Draw("hist E2");
        h2 -> Draw("same hist E2");

        legend = new TLegend(0.1,0.8,0.25,0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h1, "Signal");
        legend->AddEntry(h2, "Background");
        legend-> Draw("same"); 
        c1 -> Write(h1_names[k]);
    }



    ///SIGNAL OVER BACKROUND STUFF

    h1 = (TH1D*) f1 -> Get("4eEventLevel/hEv_HReco_M");
    h2 = (TH1D*) f2 -> Get("4eEventLevel/hEv_HReco_M");   
    

    for(Int_t i = 1; i <= h1 -> GetNbinsX(); i++){
        //std::cout << h1->GetBinLowEdge(i) << std::endl;
    }
    


    //h1 -> SetAxisRange(120, 130, "X");
    h1 -> Divide(h2);
    //h1 -> SetAxisRange(0, (h1 -> GetMaximum())*1.1, "Y");
    h1 -> SetAxisRange(0, 100, "Y");
    h1 -> Draw("hist E2");
    c1 -> Print("TestCanvas.png");


    OutputFile -> Close();
}