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

    TH1D *h1, *h2, *h_h, *h_Z, *h_Zs, *hs, *hb;
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
    h_h -> GetXaxis()->SetRangeUser(0, 130); //sets the visual range, like my old histogram goes from 0 to 150 but i wanna show till 130
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
    h_h -> GetXaxis()->SetRangeUser(0, 130); //sets the visual range, like my old histogram goes from 0 to 150 but i wanna show till 130
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

 


    //TH1 loop
    for(Int_t k = 0; k < h1_names.size() ; k++){
        h1 = (TH1D*) f1 -> Get(h1_names[k]);
        h2 = (TH1D*) f2 -> Get(h1_names[k]);    
        
        histmax = h1 -> GetMaximum();
        if(histmax < h2 -> GetMaximum()) histmax = h2 -> GetMaximum();

        h1 -> SetAxisRange(0, histmax*1.1, "Y");
        h1 -> SetLineColor(2);
        h2 -> SetLineColor(4);

        h1 -> Draw("hist E2");
        h2 -> Draw("same hist E2");

        legend = new TLegend(0.1,0.8,0.25,0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h1, "Signal");
        legend->AddEntry(h2, "Background");
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


    ///SIGNAL OVER BACKROUND STUFF

    h1 = (TH1D*) f1 -> Get("4eEventLevel/hEv_HReco_M");
    h2 = (TH1D*) f2 -> Get("4eEventLevel/hEv_HReco_M");

    THStack *hstk = new THStack("hstk","");
    hstk -> Add(h1);
    hstk -> Add(h2);
    hstk -> Draw();
    c1 -> Write("stacktest");

    

    //for(Int_t i = 1; i <= h1 -> GetNbinsX(); i++){
     //   std::cout << h1->GetBinLowEdge(i) << std::endl;
    //}
    


    //h1 -> SetAxisRange(120, 130, "X");
    h1 -> Divide(h2);
    //h1 -> SetAxisRange(0, (h1 -> GetMaximum())*1.1, "Y");
    h1 -> SetAxisRange(0, 100, "Y");
    h1 -> Draw("hist E2");
    c1 -> Print("TestCanvas.png");


    //Cut analysis//Mine

    h1 = (TH1D*) f1 -> Get("4eEventLevel/CutsAnalysis/hEvC_Logy");
    h2 = (TH1D*) f2 -> Get("4eEventLevel/CutsAnalysis/hEvC_Logy");  

    
    
    //h1 -> SetAxisRange(120, 130, "X");

    double sob;
    for(int i = 0 ; i < h1->GetNbinsX(); i++){
        std::cout << h1 -> GetBinLowEdge(i) << " " << h1->GetBinContent(i) << " " << h2->GetBinContent(i) << std::endl;
        
        sob  =  h1->GetBinContent(i) / h2->GetBinContent(i);
        if(std::isinf(sob) || std::isnan(sob)) sob = 0;
        h1 -> SetBinContent(i, sob);

    }

    
    std::cout<<h1->GetNbinsX()<<std::endl;
    //h1 -> Divide(h2);
    h1 -> SetAxisRange(0, (h1->GetBinContent(h1->GetMaximumBin()) * 1.1), "Y");
    //h1 -> SetAxisRange(0, 500, "Y");

    
    
    h1 -> Draw("hist E2");
    c1 -> Print("TestCanvas.png");
    c1 -> Write("TestCanvas");

    
    

    std::cout<<"test3"<<std::endl;


    //Cut analysis Andy's example//
    //retrieve the cut_value vector
    std::vector<double> *tmp;
    f1 -> GetObject("cut_values", tmp);
    std::vector<double> cut_values = *tmp;
    
    TH1D* h_Scan_Signif = new TH1D("h_Scan_Signif","; Log y Cut; S/sqrt(B);",cut_values.size(),cut_values.front(),cut_values.back());
    TH1D* h_Scan_Snr = new TH1D("h_Scan_Snr","; Log y Cut; S/B;",cut_values.size(),cut_values.front(),cut_values.back());

    TString histcutname;

    for(int n = 0; n < cut_values.size(); ++n){
        histcutname = "4eEventLevel/CutsAnalysis/Example/h_varycut_cut" + std::to_string(n);
        //std::cout << histcutname << std::endl;


        hs = (TH1D*) f1 -> Get(histcutname);
        hb = (TH1D*) f2 -> Get(histcutname);

        //std::cout << hs -> GetMaximum() << std::endl;
        
        int bin_low = hs ->FindBin(120);
        int bin_high= hs ->FindBin(130);

        double integral_sig = hs -> Integral(bin_low, bin_high);
        double integral_bkgd = hb -> Integral(bin_low, bin_high); 

        double signif = integral_sig / sqrt(integral_bkgd);
        double snr = integral_sig/integral_bkgd;

        int bin_thing = h_Scan_Signif -> FindBin(cut_values.at(n));

        
        if(std::isinf(signif) || std::isnan(signif)) signif = 0;
        if(std::isinf(snr) || std::isnan(snr)) snr = 0;

        
        
        std::cout << "-----------" << std::endl;
		std::cout << "Cut value = " << cut_values.at(n) << " GeV" << std::endl;
		//std::cout << "S/sqrt(B) = " << signif << std::endl;
		std::cout << "S = " <<integral_sig << ", B = " << integral_bkgd << ", S/B = " << snr << std::endl;
        //std::cout << "is inf? " << std::isinf(snr) << std::endl;
        

        

        h_Scan_Signif -> SetBinContent(bin_thing, signif);
        h_Scan_Snr -> SetBinContent(bin_thing, snr);
        
    }
    std::cout<<h_Scan_Snr -> GetNbinsX() << std::endl;

    h_Scan_Signif->Write();
	h_Scan_Snr->Write();


    OutputFile -> Close();
}