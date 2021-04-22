void fit_bgsig()
{
        // style setting
        gROOT->SetStyle("Plain");
        gStyle->SetOptFit(111111);
        gStyle->SetFrameBorderMode(0);
        gStyle->SetFillColor(0);
        gStyle->SetCanvasColor(0);
        gStyle->SetFrameFillColor(0);
        gStyle->SetCanvasBorderMode(0);

        // create a random number generator
        gRandom = new TRandom3();

        // create a histogram 
        TH1D * hist = new TH1D("data", ";x;N", 20, 0.0, 100.0);

        // fill in signal
        for (int i = 0; i < 100; ++i)
                hist->Fill(gRandom->Gaus(65.0, 5.0));
        // fill background
        for (int i = 0; i < 100; ++i)
                hist->Fill(100.0 * gRandom->Uniform());

        // define a fit function = gauss + pol0 this is the important line for the histogram fitting
        TF1 * f1 = new TF1("myfunc", "[0] / sqrt(2.0 * TMath::Pi()) / [2] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3]", 0, 100);

        //set parameter start values (mandatory). 
        f1->SetParNames("Constant","Mean","Sigma","Flat");
        f1->SetParameters(700.,hist->GetMean(),hist->GetRMS(),2.);
        f1->SetParLimits(0,  100.0, 700.0);
        f1->SetParLimits(1, 50.0,  90.0);
        f1->SetParLimits(2,  0.1,  10.0);
        f1->SetParLimits(3,  0.0,   10.0);
        f1->SetLineWidth(2);
        f1->SetLineColor(2);

        // create a canvas to draw the histogram
        TCanvas * c1= new TCanvas("c1", "fit sig+bg",5,5,800,600);
        // perform fit
        hist->Fit("myfunc");
        hist->Draw();
        hist->SaveAs("fit_bgsig.eps");

}