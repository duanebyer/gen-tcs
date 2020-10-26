#include <iostream>
#include <string>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TTree.h>

// Convert from pb to cm^2, and from s^-1 to hr^-1.
double const CROSS_SECTION_TO_HOURLY_RATE = 1e-36 * 60 * 60;

// Convenience program for integrating over the background between two values.
int main(int argc, char** argv) {
	if (argc != 5) {
		std::cout << "Usage: integrate <ROOT file> <lum.> <M_ll lower> <M_ll upper>" << std::endl;
		return 1;
	}

	TFile* file_in = new TFile(argv[1], "READ");
	if (file_in->IsZombie()) {
		std::cerr << "Error opening input file " << argv[1] << std::endl;
		return 1;
	}
	Double_t luminosity = std::stod(argv[2]);
	Double_t M_ll_lower = std::stod(argv[3]);
	Double_t M_ll_upper = std::stod(argv[4]);

	TTree* events = (TTree*) file_in->Get("tr1");
	TParameter<Int_t>* param_Nsim = (TParameter<Int_t>*) file_in->Get("Nsim");
	TParameter<Int_t>* param_Nsim_tot = (TParameter<Int_t>*) file_in->Get("Nsim_tot");
	Int_t Nsim = param_Nsim->GetVal();
	Int_t Nsim_tot = param_Nsim_tot->GetVal();

	TLorentzVector* L_em = nullptr;
	TLorentzVector* L_ep = nullptr;
	Double_t psf;
	Double_t flux_factor;
	Double_t flux_brem;
	Double_t crs_BH;
	events->SetBranchAddress("L_em", &L_em);
	events->SetBranchAddress("L_ep", &L_ep);
	events->SetBranchAddress("psf", &psf);
	events->SetBranchAddress("flux_factor", &flux_factor);
	events->SetBranchAddress("flux_brem", &flux_brem);
	events->SetBranchAddress("crs_BH", &crs_BH);

	Int_t event_count = events->GetEntries();
	Double_t total_rate = 0.;
	for (Int_t event_idx = 0; event_idx < event_count; ++event_idx) {
		events->GetEntry(event_idx);
		Double_t M_ll = (*L_em + *L_ep).M();
		if (M_ll >= M_ll_lower && M_ll <= M_ll_upper) {
			Double_t rate = crs_BH * psf * (flux_factor + flux_brem) * luminosity * CROSS_SECTION_TO_HOURLY_RATE / Nsim_tot;
			total_rate += rate;
		}
	}

	std::cout << "Total rate: " << total_rate << " counts / hour" << std::endl;
	return 0;
}

