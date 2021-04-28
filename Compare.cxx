#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
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
    TString filename3 = "bakgndZ.root";


    TFile *f1 = TFile::Open(filename1);
    TFile *f2 = TFile::Open(filename2);
    TFile *f3 = TFile::Open(filename3);

    TFile *OutputFile = new TFile("Compare.root", "recreate");
    OutputFile -> cd();
    //TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,1680,1050); //original
    TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,640,480); //for quicklook
    gStyle -> SetOptStat(0);

    TH1D *h1, *h2, *h3, *herr, *h_h, *h_Z, *h_Zs, *hs, *hb;
    TLegend *legend;
    Double_t histmax;
    Float_t rightmax, scale;
    TGaxis *axis;

    
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

    h1_names.push_back("4eEventLevel/CutsAnalysis/hEvC_Zstar");
    h1_names.push_back("4eEventLevel/CutsAnalysis/hEvC_Logy");


    //Combining Z, Z*, and m4l
    h_h  = (TH1D*) f1 -> Get("4eEventLevel/hEv_HReco_M");
    h_Z  = (TH1D*) f1 -> Get("4eEventLevel/hEv_ZReco_M");
    h_Zs = (TH1D*) f1 -> Get("4eEventLevel/hEv_ZstarReco_M");

    histmax = h_h -> GetMaximum();
    if(histmax < h_Z -> GetMaximum()) histmax = h_Z -> GetMaximum();
    if(histmax < h_Zs -> GetMaximum()) histmax = h_Zs -> GetMaximum();

    //note that here we're applying the formatting on the _first_ histogram we wanna stick on tcanvas
    h_h -> SetAxisRange(0, histmax*1.1, "Y");
    h_h -> SetTitle("Reconstruction of Higgs, Z, and Z* Bosons"); //sets title of the histogram ur gonna print
    h_h ->GetXaxis()->SetTitle("Mass (GeV)");  //sets x axis title
    h_h ->GetYaxis()->SetTitle("Number of Events"); //sets y axis title
    h_h -> GetXaxis()->SetRangeUser(0, 200); //sets the visual range, like my old histogram goes from 0 to 150 but i wanna show till 130
    h_h -> SetLineColor(kRed); //changes the bars of the histogram so its red
    h_h -> SetStats(kFALSE); //removes the annoying 'mean median and entries' box thing
    h_h -> Draw("hist E2"); //Draws the histogram using the option hist and E2
    //The Draw command takes h_h with its applied options and sticks it into the TCanvas canvas

    h_Z -> SetLineColor(kBlue); //this is my second histogram
    h_Z -> Draw("hist same E2"); //NOTE THE 'SAME' BIT, draws the histogram on the same canvas

    h_Zs->SetLineColor(kGreen); //third histogram
    h_Zs -> Draw("same hist E2");

    legend = new TLegend(0.1,0.7,0.35,0.9); // Initialize the TLegend object that you place within the TCanvas object
    //numbers are (x1, y1, x2, y2), bottom left and top right corner (i think)
    legend->SetHeader("Particles","C"); // option "C" allows to center the header. 'Particles' is the Title of the legend
    legend->AddEntry(h_h, "Higgs"); //Adds the entry and the label
    legend->AddEntry(h_Z,"Z Boson");
    legend->AddEntry(h_Zs,"Z* Boson");
    legend-> Draw("same");  //Draws on same histogram

    c1 -> Write("ReconstructedMassSignal");

    h_h  = (TH1D*) f2 -> Get("4eEventLevel/hEv_HReco_M");
    h_Z  = (TH1D*) f2 -> Get("4eEventLevel/hEv_ZReco_M");
    h_Zs = (TH1D*) f2 -> Get("4eEventLevel/hEv_ZstarReco_M");

    histmax = h_h -> GetMaximum();
    if(histmax < h_Z -> GetMaximum()) histmax = h_Z -> GetMaximum();
    if(histmax < h_Zs -> GetMaximum()) histmax = h_Zs -> GetMaximum();

    //note that here we're applying the formatting on the _first_ histogram we wanna stick on tcanvas
    h_h -> SetAxisRange(0, histmax*1.1, "Y");
    h_h -> SetTitle("Reconstruction of Higgs, Z, and Z* Bosons"); //sets title of the histogram ur gonna print
    h_h ->GetXaxis()->SetTitle("Mass (GeV)");  //sets x axis title
    h_h ->GetYaxis()->SetTitle("Number of Events"); //sets y axis title
    h_h -> GetXaxis()->SetRangeUser(0, 200); //sets the visual range, like my old histogram goes from 0 to 150 but i wanna show till 130
    h_h -> SetLineColor(kRed); //changes the bars of the histogram so its red
    h_h -> SetStats(kFALSE); //removes the annoying 'mean median and entries' box thing
    h_h -> Draw("hist E2"); //Draws the histogram using the option hist and E2
    //The Draw command takes h_h with its applied options and sticks it into the TCanvas canvas

    h_Z -> SetLineColor(kBlue); //this is my second histogram
    h_Z -> Draw("hist same E2"); //NOTE THE 'SAME' BIT, draws the histogram on the same canvas

    h_Zs->SetLineColor(kGreen); //third histogram
    h_Zs -> Draw("same hist E2");

    legend = new TLegend(0.1,0.8,0.25,0.9); // Initialize the TLegend object that you place within the TCanvas object
    //numbers are (x1, y1, x2, y2), bottom left and top right corner (i think)
    legend->SetHeader("Particles","C"); // option "C" allows to center the header. 'Particles' is the Title of the legend
    legend->AddEntry(h_h, "Higgs"); //Adds the entry and the label
    legend->AddEntry(h_Z,"Z Boson");
    legend->AddEntry(h_Zs,"Z* Boson");
    legend-> Draw("same");  //Draws on same histogram

    c1 -> Write("ReconstructedMassBakgnd");

    h_h  = (TH1D*) f3 -> Get("4eEventLevel/hEv_HReco_M");
    h_Z  = (TH1D*) f3 -> Get("4eEventLevel/hEv_ZReco_M");
    h_Zs = (TH1D*) f3 -> Get("4eEventLevel/hEv_ZstarReco_M");

    histmax = h_h -> GetMaximum();
    if(histmax < h_Z -> GetMaximum()) histmax = h_Z -> GetMaximum();
    if(histmax < h_Zs -> GetMaximum()) histmax = h_Zs -> GetMaximum();

    //note that here we're applying the formatting on the _first_ histogram we wanna stick on tcanvas
    h_h -> SetAxisRange(0, histmax*1.1, "Y");
    h_h -> SetTitle("Reconstruction of Higgs, Z, and Z* Bosons"); //sets title of the histogram ur gonna print
    h_h ->GetXaxis()->SetTitle("Mass (GeV)");  //sets x axis title
    h_h ->GetYaxis()->SetTitle("Number of Events"); //sets y axis title
    h_h -> GetXaxis()->SetRangeUser(0, 200); //sets the visual range, like my old histogram goes from 0 to 150 but i wanna show till 130
    h_h -> SetLineColor(kRed); //changes the bars of the histogram so its red
    h_h -> SetStats(kFALSE); //removes the annoying 'mean median and entries' box thing
    h_h -> Draw("hist E2"); //Draws the histogram using the option hist and E2
    //The Draw command takes h_h with its applied options and sticks it into the TCanvas canvas

    h_Z -> SetLineColor(kBlue); //this is my second histogram
    h_Z -> Draw("hist same E2"); //NOTE THE 'SAME' BIT, draws the histogram on the same canvas

    h_Zs->SetLineColor(kGreen); //third histogram
    h_Zs -> Draw("same hist E2");

    legend = new TLegend(0.1,0.8,0.25,0.9); // Initialize the TLegend object that you place within the TCanvas object
    //numbers are (x1, y1, x2, y2), bottom left and top right corner (i think)
    legend->SetHeader("Particles","C"); // option "C" allows to center the header. 'Particles' is the Title of the legend
    legend->AddEntry(h_h, "Higgs"); //Adds the entry and the label
    legend->AddEntry(h_Z,"Z Boson");
    legend->AddEntry(h_Zs,"Z* Boson");
    legend-> Draw("same");  //Draws on same histogram

    c1 -> Write("ReconstructedMassBakgndZ");

 
    h_h  = (TH1D*) f1 -> Get("4eEventLevel/hEv_H_eta");
    h_Z  = (TH1D*) f1 -> Get("4eEventLevel/hEv_jet_eta");
    h_Zs = (TH1D*) f1 -> Get("4eEventLevel/hEv_nue_eta_wicuts");

    histmax = h_h -> GetMaximum();
    if(histmax < h_Z -> GetMaximum()) histmax = h_Z -> GetMaximum();
    if(histmax < h_Zs -> GetMaximum()) histmax = h_Zs -> GetMaximum();

    //note that here we're applying the formatting on the _first_ histogram we wanna stick on tcanvas
    h_h -> SetAxisRange(0, histmax*1.1, "Y");
    h_h -> SetTitle("Eta Distribution of Higgs, Scattered Quark (jet) and Scattered Lepton (neutrino)"); //sets title of the histogram ur gonna print
    h_h ->GetXaxis()->SetTitle("Eta");  //sets x axis title
    h_h ->GetYaxis()->SetTitle("Number of Events"); //sets y axis title
    h_h -> SetLineColor(kRed); //changes the bars of the histogram so its red
    h_h -> SetStats(kFALSE); //removes the annoying 'mean median and entries' box thing
    h_h -> Draw("hist E2"); //Draws the histogram using the option hist and E2

    h_Z -> SetLineColor(kGreen); //this is my second histogram
    h_Z -> Draw("hist same E2"); //NOTE THE 'SAME' BIT, draws the histogram on the same canvas

    h_Zs->SetLineColor(kBlue); //third histogram
    h_Zs -> Draw("same hist E2");

    legend = new TLegend(0.1,0.7,0.35,0.9); // Initialize the TLegend object that you place within the TCanvas object
    //numbers are (x1, y1, x2, y2), bottom left and top right corner (i think)
    legend->SetHeader("Particles","C"); // option "C" allows to center the header. 'Particles' is the Title of the legend
    legend->AddEntry(h_h, "Higgs"); //Adds the entry and the label
    legend->AddEntry(h_Z,"Jet");
    legend->AddEntry(h_Zs,"Neutrino");
    legend-> Draw("same");  //Draws on same histogram

    c1 -> Write("EtaPlot");

    //TH1 loop
    for(Int_t k = 0; k < h1_names.size() ; k++){
        h1 = (TH1D*) f1 -> Get(h1_names[k]);
        h2 = (TH1D*) f2 -> Get(h1_names[k]);
        h3 = (TH1D*) f3 -> Get(h1_names[k]);    
        
        histmax = h1 -> GetMaximum();
        if(histmax < h2 -> GetMaximum()) histmax = h2 -> GetMaximum();
        if(histmax < h3 -> GetMaximum()) histmax = h3 -> GetMaximum();

        h1 -> SetAxisRange(0, histmax*1.1, "Y");
        h1 -> SetLineColor(2);
        h2 -> SetLineColor(4);
        h3 -> SetLineColor(8);

        h1 -> Draw("hist E2");
        h2 -> Draw("same hist E2");
        h3 -> Draw("same hist E2");

        legend = new TLegend(0.1,0.8,0.25,0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h1, "Signal");
        legend->AddEntry(h2, "ZZ Background");
        legend->AddEntry(h3, "Z Background");
        legend-> Draw("same"); 
        c1 -> Write(h1_names[k] + "_actual");

        
        h1->Draw("hist E2");
        c1 -> Update();
        rightmax = 1.1*h2->GetBinContent(h2->GetMaximumBin());
        scale = gPad->GetUymax()/rightmax;
        //std::cout<<rightmax<<std::endl;
        h2->Scale(scale); //SCALE
        axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), 
                          gPad->GetUxmax(), gPad->GetUymax(),
                          0, rightmax,  510 , "+L");
        axis -> SetLineColor(4);
        axis -> SetLabelColor(4);
        h2->Draw("same hist E2");
        axis->Draw();
        c1 -> Write(h1_names[k] + "_scale");
        h2 -> Scale(1/scale);
    
    }


    ///Stacked Histogram Stuff
    std::vector<TString> h1St_names;
    h1St_names.push_back("4eEventLevel/Smearing/hEv_HReco_M");
    h1St_names.push_back("4eEventLevel/Smearing/hEv_HReco_M_S");

    


    


    


    //h1 -> SetAxisRange(120, 130, "X");
    h1 -> Divide(h2);
    //h1 -> SetAxisRange(0, (h1 -> GetMaximum())*1.1, "Y");
    h1 -> SetAxisRange(0, 100, "Y");
    h1 -> Draw("hist E2");
    //c1 -> Print("TestCanvas.png");


    //Cut analysis//Mine

    h1 = (TH1D*) f1 -> Get("4eEventLevel/CutsAnalysis/hEvC_Logy");
    h2 = (TH1D*) f2 -> Get("4eEventLevel/CutsAnalysis/hEvC_Logy");  

    

    
    
    //h1 -> SetAxisRange(120, 130, "X");

    double sob;
    for(int i = 0 ; i < h1->GetNbinsX(); i++){
        //std::cout << h1 -> GetBinLowEdge(i) << " " << h1->GetBinContent(i) << " " << h2->GetBinContent(i) << std::endl;
        
        sob  =  h1->GetBinContent(i) / h2->GetBinContent(i);
        if(std::isinf(sob) || std::isnan(sob)) sob = 0;
        h1 -> SetBinContent(i, sob);

    }

    
    //std::cout<<h1->GetNbinsX()<<std::endl;
    //h1 -> Divide(h2);
    h1 -> SetAxisRange(0, (h1->GetBinContent(h1->GetMaximumBin()) * 1.1), "Y");
    //h1 -> SetAxisRange(0, 500, "Y");

    
    
    h1 -> Draw("hist E2");
    c1 -> Write("TestCanvas");





    //Smearing Framework
    std::vector<TString> h1S_names;
    h1S_names.push_back("4eEventLevel/Smearing/hEvS_e_pt");
    h1S_names.push_back("4eEventLevel/Smearing/hEvS_e_E");
    h1S_names.push_back("4eEventLevel/Smearing/hEv_HReco_M");

    for(Int_t k = 0; k < h1S_names.size(); k++){
        h1->Reset();
        h2->Reset();
        h1 = (TH1D*) f1 -> Get(h1S_names[k]);
        h2 = (TH1D*) f1 -> Get(h1S_names[k] + "_S");
        histmax = h1 -> GetMaximum();
        if(histmax < h2 -> GetMaximum()) histmax = h2 -> GetMaximum();

        h1 -> SetAxisRange(0, histmax*1.1, "Y");
        h1 -> SetLineColor(2);
        h1 -> SetFillColor(2);
        h2 -> SetLineColor(4);

        //h2-> SetFillColor(7);
        //h2-> SetFillColorAlpha(7,0);


        h1 -> Draw("hist");
        h2 -> Draw("same hist");

        legend = new TLegend(0.1,0.8,0.25,0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h1, "Signal");
        legend->AddEntry(h2, "Smearing");
        legend-> Draw("same"); 
        c1 -> Write(h1S_names[k]);
    }

    for(Int_t k = 0; k < h1St_names.size(); k++){

        h1 = (TH1D*) f1 -> Get(h1St_names[k]);
        h2 = (TH1D*) f2 -> Get(h1St_names[k]);
        h3 = (TH1D*) f3 -> Get(h1St_names[k]);

        h1 -> SetFillColor(38);
        h2 -> SetFillColor(2);
        h3 -> SetFillColor(46);

        h1 -> SetLineColor(38);
        h2 -> SetLineColor(2);
        h3 -> SetLineColor(46);
        
        THStack *hstk = new THStack("hstk","");
        hstk -> Add(h2);
        hstk -> Add(h3);
        hstk -> Add(h1);
        hstk -> Draw("HIST");

        legend = new TLegend(0.1,0.7,0.35,0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h1, "Signal");
        legend->AddEntry(h2, "ZZ Background");
        legend->AddEntry(h3, "Z Background");
        legend-> Draw("same"); 

        hstk -> SetTitle("Stacked Histogram of Mass Distribution of H, ZZ*, and Z/Z*"); //sets title of the histogram ur gonna print
        hstk ->GetXaxis()->SetTitle("M_{4l} (GeV)");  //sets x axis title
        hstk ->GetYaxis()->SetTitle("Events"); //sets y axis title
        TString titlestack = Form("%d",k);
        c1 -> Write("stacktest" + titlestack);

        




        //Finding out number of events
        TAxis *xaxisstk = hstk->GetXaxis();
        Int_t binxstklow = xaxisstk->FindBin(100);
        Int_t binxstkhigh = xaxisstk->FindBin(150);
        Int_t binxstk = xaxisstk->FindBin(125);
        TH1 *last = (TH1*)hstk->GetStack()->Last(); 

        Int_t lengthstack = hstk->GetStack()->GetLast();
        

        //std::cout<< lengthstack << std::endl;
        std::cout << last->Integral(binxstklow, binxstkhigh) << std::endl; //total
        std::cout << h1->Integral(binxstklow, binxstkhigh) << std::endl; //higgs
        std::cout << h2->Integral(binxstklow, binxstkhigh) << std::endl; //ZZ
        std::cout << h3->Integral(binxstklow, binxstkhigh) << std::endl; //Z
        std::cout << " " << std::endl;

        
        //new thing needs to be in steps of 3
        Int_t binstart = xaxisstk->FindBin(125);
        Int_t binvarhigh, binvarlow;
        float totalvar, higgsvar, zzvar, zvar; 

        float testnumber;
        
        if(k ==1 ){        
            for(Int_t k = 0; k <= 75; k = k + 3){
                binvarhigh = xaxisstk->FindBin(125 + k);
                binvarlow  = xaxisstk->FindBin(125 - k);
                //i now have all the information i need
                totalvar = last->Integral(binvarlow, binvarhigh); 
                higgsvar = h1  ->Integral(binvarlow, binvarhigh);
                zzvar    = h2  ->Integral(binvarlow, binvarhigh);
                zvar     = h3  ->Integral(binvarlow, binvarhigh);

                testnumber = TMath::Sqrt(totalvar +zzvar + zvar) / (totalvar - zzvar - zvar);//precision we can measure 
                std::cout <<
                    testnumber*100 << "%    " << 
                    last->GetBinCenter(binvarlow) << " " << 
                    last->GetBinCenter(binvarhigh) << " "<< std::endl;
                
            }

            
        }
        
    }

    
    
    


    OutputFile -> Close();
}