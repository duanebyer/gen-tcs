#include <fstream>
#include <iostream>
#include <string>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TTree.h>

// Converts a file from a .root format to a .dat format.
int main(int argc, char** argv) {
	if (argc != 3) {
		std::cout << "Usage: root-to-dat <file in> <file out>" << std::endl;
		return 1;
	}
	TFile* file_in = new TFile(argv[1], "READ");
	if (file_in->IsZombie()) {
		std::cerr << "Error opening file to read " << argv[1] << std::endl;
		return 1;
	}
	std::ofstream file_out(argv[2]);
	if (!file_out.good()) {
		std::cerr << "Error opening file to write " << argv[2] << std::endl;
		return 1;
	}

	TTree* events = (TTree*) file_in->Get("tr1");
	TParameter<Double_t>* param_Eb = (TParameter<Double_t>*) file_in->Get("Eb");
	TParameter<Int_t>* param_Nsim = (TParameter<Int_t>*) file_in->Get("Nsim");
	TParameter<Int_t>* param_Nsim_tot = (TParameter<Int_t>*) file_in->Get("Nsim_tot");

	Double_t Eb = param_Eb->GetVal();
	Int_t Nsim = param_Nsim->GetVal();
	Int_t Nsim_tot = param_Nsim_tot->GetVal();
	file_out << "Ebeam=" << Eb << "GeV" << ",";
	file_out << "Nsim= " << Nsim << std::endl;

	TLorentzVector* L_em = nullptr;
	TLorentzVector* L_ep = nullptr;
	TLorentzVector* L_prot = nullptr;
	TLorentzVector* L_g = nullptr;
	Double_t Eg;
	Double_t psf;
	Double_t flux_factor;
	Double_t flux_brem;
	Double_t crs_BH;
	events->SetBranchAddress("L_em", &L_em);
	events->SetBranchAddress("L_ep", &L_ep);
	events->SetBranchAddress("L_prot", &L_prot);
	events->SetBranchAddress("L_g", &L_g);
	events->SetBranchAddress("Eg", &Eg);
	events->SetBranchAddress("psf", &psf);
	events->SetBranchAddress("flux_factor", &flux_factor);
	events->SetBranchAddress("flux_brem", &flux_brem);
	events->SetBranchAddress("crs_BH", &crs_BH);

	Int_t event_count = events->GetEntries();
	for (Int_t event_idx = 0; event_idx < event_count; ++event_idx) {
		events->GetEntry(event_idx);
		Double_t N_ratio = (Double_t) Nsim / (Double_t) Nsim_tot;
		// Convert the weight to use GeV^-2 instead of pb.
		Double_t weight = crs_BH * psf * (flux_brem + flux_factor) * N_ratio * 2.56819e-9;
		file_out << weight << std::endl;
		file_out << "q:\t" << L_g->X() << '\t' << L_g->Y() << '\t' << L_g->Z() << '\t' << L_g->T() << std::endl;
		file_out << "e+:\t" << L_ep->X() << '\t' << L_ep->Y() << '\t' << L_ep->Z() << '\t' << L_ep->T() << std::endl;
		file_out << "e-:\t" << L_em->X() << '\t' << L_em->Y() << '\t' << L_em->Z() << '\t' << L_em->T() << std::endl;
		file_out << "p:\t" << L_prot->X() << '\t' << L_prot->Y() << '\t' << L_prot->Z() << '\t' << L_prot->T() << std::endl;
	}

	return 0;
}

