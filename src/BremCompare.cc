#include <iostream>

#include <KinFunctions.h>

#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>

Double_t epa_flux(Double_t* x, Double_t* p) {
	Double_t Eb = p[0];
	Double_t q2_cut = p[1];
	return KinFuncs::N_EPA(Eb, x[0], q2_cut);
}

Double_t brem_flux(Double_t* x, Double_t* p) {
	Double_t Eb = p[0];
	Double_t l_targ = p[1];
	return KinFuncs::Brem_Approx(x[0], Eb, l_targ);
}

int main(int argc, char** argv) {
	if (argc != 4) {
		std::cout << "Usage: brem-compare <Eb> <Q2 cut> <targ. length>" << std::endl;
		return 1;
	}
	Double_t Eb = std::stod(argv[1]);
	Double_t q2_cut = std::stod(argv[2]);
	Double_t l_targ = std::stod(argv[3]);

	TApplication* app = new TApplication("app", 0, 0);
	TF1* epa = new TF1("N_EPA", &epa_flux, 0, Eb, 2);
	TF1* brem = new TF1("brem", &brem_flux, 0, Eb, 2);
	epa->SetParameter(0, Eb);
	epa->SetParameter(1, q2_cut);
	brem->SetParameter(0, Eb);
	brem->SetParameter(1, l_targ);
	epa->SetLineColor(1);
	epa->GetXaxis()->SetTitle("virtual photon energy (GeV)");
	epa->GetYaxis()->SetTitle("flux factor");
	brem->SetLineColor(2);

	TCanvas* canvas = new TCanvas();
	canvas->SetLogy();
	TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
	legend->AddEntry(epa, "epa flux factor", "l");
	legend->AddEntry(brem, "brem flux factor", "l");
	epa->Draw();
	brem->Draw("same");
	legend->Draw();

	app->Run();
}
