#include <fstream>
#include <iostream>
#include <string>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

// Converts a file from a .root format to a .dat format.
int main(int argc, char** argv) {
	if (argc != 3) {
		std::cout << "Usage: RootToDat <file in> <file out>" << std::endl;
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

	TLorentzVector L_em;
	TLorentzVector L_ep;
	TLorentzVector L_prot;
	Double_t E_g;
	Double_t psf;
	Double_t flux_factor;
	Double_t flux_brem;
	Double_t crs_BH;
	events->SetBranchAddress("L_em", &L_em);
	events->SetBranchAddress("L_ep", &L_ep);
	events->SetBranchAddress("L_prot", &L_prot);
	events->SetBranchAddress("E_g", &E_g);
	events->SetBranchAddress("psf", &psf);
	events->SetBranchAddress("flux_factor", &flux_factor);
	events->SetBranchAddress("flux_brem", &flux_brem);
	events->SetBranchAddress("crs_BH", &crs_BH);

	Int_t event_count = events->GetEntries();
	for (Int_t event_idx = 0; event_idx < event_count; ++event_idx) {
		events->GetEntry(event_idx);
		Double_t N_ratio = (Double_t) Nsim / (Double_t) Nsim_tot;
		Double_t weight = crs_BH * psf * (flux_brem + flux_factor) * N_ratio;
		TLorentzVector L_q(E_g, 0., 0., E_g);
		file_out << weight << std::endl;
		file_out << "q:\t" << L_q.X() << '\t' << L_q.Y() << '\t' << L_q.Z() << '\t' << L_q.T() << std::endl;
		file_out << "e+:\t" << L_ep.X() << '\t' << L_ep.Y() << '\t' << L_ep.Z() << '\t' << L_ep.T() << std::endl;
		file_out << "e-:\t" << L_em.X() << '\t' << L_em.Y() << '\t' << L_em.Z() << '\t' << L_em.T() << std::endl;
		file_out << "p:\t" << L_prot.X() << '\t' << L_prot.Y() << '\t' << L_prot.Z() << '\t' << L_prot.T() << std::endl;
	}

	return 0;
}

