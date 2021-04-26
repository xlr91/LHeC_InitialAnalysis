void f12() {
   Double_t angle = TMath::Pi()/4.;
   Double_t cosphi = TMath::Cos(angle);
   Double_t sinphi = TMath::Sin(angle);
   
   TH2F *h2 = new TH2F("h2","h2",100,0,0.8,100,0,0.8);
   TRandom r;     
   for (Int_t i=0;i<1000;i++) {
      Double_t x = r.Uniform(0,1);
      Double_t y = r.Gaus(0,0.02);
      Double_t u = cosphi*x -sinphi*y;
      Double_t v = sinphi*x +cosphi*y;
      h2->Fill(u,v);
   }
   TF1 *f1 = new TF1("f1","[0]+[1]*x",0,1);
   f1->SetParameters(0.,1.);
   f1->SetLineColor(kRed);
   h2->Fit(f1);
   h2->Draw();
   f1->Draw("same");
}