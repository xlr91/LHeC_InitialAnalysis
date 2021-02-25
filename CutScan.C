
void CutScan() {

	TFile* file = new TFile("CutOptTest.root","recreate");

	std::vector<TH1D*> h_varycut_signal;
	std::vector<TH1D*> h_varycut_bkgd;

	int nCuts = 10;
	double min_cut = 10.0; // GeV
	double step_cut = 2.0; // GeV

	std::vector<double> cut_values;

	for(int n = 0; n < nCuts; ++n) {

		TString suffix = "_cut";
		suffix += n;

		TH1D* h_new_signal = new TH1D("h_varycut_signal"+suffix,"; m_{4l} [GeV]; Events / 4 GeV",50,50,250);
		TH1D* h_new_bkgd = new TH1D("h_varycut_bkgd"+suffix,"; m_{4l} [GeV]; Events / 4 GeV",50,50,250);

		h_varycut_signal.push_back(h_new_signal);
		h_varycut_bkgd.push_back(h_new_bkgd);

		cut_values.push_back( min_cut + n*step_cut );

	}

	std::cout << h_varycut_signal.size() << " " << h_varycut_bkgd.size() << std::endl;

	TRandom3* m_rand = new TRandom3(4358);

	const int nEvents_Signal = 10000;
	const int nEvents_Bkgd = 10000;

	// Loop over signal events
	for(int i = 0; i < nEvents_Signal; ++i) {

		const double lepton_pT = 10.0 + 2.0*m_rand->Uniform(0.0,30.0);
		const double m_4l = m_rand->Gaus(125,5.0);

		// Loop over cuts to check if this event passes the all
		for(int n = 0; n < nCuts; ++n) {

			if( lepton_pT > cut_values.at(n) ) {
				h_varycut_signal.at(n)->Fill(m_4l);
			}

		}

	} // Loop over signal events

	// Loop over bkgd. events
	for(int i = 0; i < nEvents_Bkgd; ++i) {

		const double lepton_pT = 10.0 + 1.5*m_rand->Uniform(0.0,30.0);
		const double m_4l = m_rand->Uniform(50,250);

		// Loop over cuts to check if this event passes the all
		for(int n = 0; n < nCuts; ++n) {

			if( lepton_pT > cut_values.at(n) ) {
				h_varycut_bkgd.at(n)->Fill(m_4l);
			}

		}

	} // Loop over signal events

	TH1D* h_Scan_Signif = new TH1D("h_Scan_Signif","; Lepton pT Cut [GeV]; S/sqrt(B);",nCuts,min_cut,min_cut + nCuts*step_cut);
	TH1D* h_Scan_Snr = new TH1D("h_Scan_Snr","; Lepton pT Cut [GeV]; S/B;",nCuts,min_cut,min_cut + nCuts*step_cut);

	// Loop over cut values
	for(int n = 0; n < nCuts; ++n) {

		int bin_low = h_varycut_signal.at(n)->FindBin(110);
		int bin_high = h_varycut_signal.at(n)->FindBin(140);

		double integral_sig = h_varycut_signal.at(n)->Integral(bin_low,bin_high);
		double integral_bkgd = h_varycut_bkgd.at(n)->Integral(bin_low,bin_high);

		double signif = integral_sig/sqrt(integral_bkgd);
		double snr = integral_sig/integral_bkgd;

		std::cout << "-----------" << std::endl;
		std::cout << "Cut value = " << cut_values.at(n) << " GeV" << std::endl;
		std::cout << "S/sqrt(B) = " << signif << std::endl;
		std::cout << "S/B = " << snr << std::endl;

		int bin_pt = h_Scan_Signif->FindBin(cut_values.at(n));

		h_Scan_Signif->SetBinContent(bin_pt,signif);
		h_Scan_Snr->SetBinContent(bin_pt,snr);

	}

	file->cd();

	h_Scan_Signif->Write();
	h_Scan_Snr->Write();

	for(int n = 0; n < nCuts; ++n) {

		h_varycut_signal.at(n)->Write();
		h_varycut_bkgd.at(n)->Write();

	}

	file->Close();

}