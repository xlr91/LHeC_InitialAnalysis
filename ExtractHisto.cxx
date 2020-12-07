#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <TCanvas.h>
#include "TF1.h"
#include "TEfficiency.h"


void ExtractHisto(){
    TString filename = "output.root";
    TFile *f = TFile::Open(filename);
    TCanvas *c1 = new TCanvas("c1"," Graphs ",50,50,1680,1050);
    TString current_histogram;
    TH1D *h1;
    Double_t ymax;
    Double_t xval;
    Double_t etalim = -5.0;


    current_histogram = "h_e_eta";
    h1 = (TH1D*) f -> Get(current_histogram);

    xval = etalim;
    ymax = h1 -> GetMaximum();
    TLine *line= new TLine(xval,0,xval,ymax);
    line -> SetLineColor(kRed);

    h1 -> Draw();
    line -> Draw("same");
    c1 -> Print("Graphs/" + current_histogram + ".png");


    current_histogram = "h_Jet_eta";
    h1 = (TH1D*) f -> Get(current_histogram);

    xval = etalim;
    ymax = h1 -> GetMaximum();

    
    TLine *line1= new TLine(xval,0,xval,ymax);
    line1 -> SetLineColor(kRed);
    
    h1 -> Draw();
    line1 -> Draw("same");
    c1 -> Print("Graphs/" + current_histogram + ".png");









}