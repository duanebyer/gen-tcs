#include <iostream>

#include <TFile.h>
#include <TParameter.h>
#include <TTree.h>

// Utility for showing the simulation parameters used to generate a ROOT file.
int main(int argc, char** argv) {
	if (argc != 2) {
		std::cerr << "Usage: params-inspect <ROOT file>" << std::endl;
		return 1;
	}
	TFile* file_in = new TFile(argv[1], "READ");
	if (file_in->IsZombie()) {
		std::cerr << "Error opening input file " << argv[1] << std::endl;
		return 1;
	}

	TParameter<Int_t>* param_seed = (TParameter<Int_t>*) file_in->Get("seed");
	TParameter<Int_t>* param_Nsim = (TParameter<Int_t>*) file_in->Get("Nsim");
	TParameter<Int_t>* param_Nsim_tot = (TParameter<Int_t>*) file_in->Get("Nsim_tot");
	TParameter<Int_t>* param_target = (TParameter<Int_t>*) file_in->Get("target");
	TParameter<Double_t>* param_Eb = (TParameter<Double_t>*) file_in->Get("Eb");
	TParameter<Double_t>* param_tMin = (TParameter<Double_t>*) file_in->Get("tMin");
	TParameter<Double_t>* param_tMax = (TParameter<Double_t>*) file_in->Get("tMax");
	TParameter<Double_t>* param_EgMin = (TParameter<Double_t>*) file_in->Get("EgMin");
	TParameter<Double_t>* param_EgMax = (TParameter<Double_t>*) file_in->Get("EgMax");
	TParameter<Double_t>* param_MinvMax = (TParameter<Double_t>*) file_in->Get("MinvMax");
	TParameter<Double_t>* param_MinvMin = (TParameter<Double_t>*) file_in->Get("MinvMin");
	TParameter<Double_t>* param_q2Cut = (TParameter<Double_t>*) file_in->Get("q2Cut");
	TParameter<Double_t>* param_target_l = (TParameter<Double_t>*) file_in->Get("target_l");
	TParameter<Double_t>* param_target_off = (TParameter<Double_t>*) file_in->Get("target_off");
	TParameter<Double_t>* param_Dterm = (TParameter<Double_t>*) file_in->Get("Dterm");
	TParameter<Double_t>* param_Theta_elec_min = (TParameter<Double_t>*) file_in->Get("Theta_elec_min");
	TParameter<Double_t>* param_Theta_elec_max = (TParameter<Double_t>*) file_in->Get("Theta_elec_max");
	TParameter<Double_t>* param_P_elec_min = (TParameter<Double_t>*) file_in->Get("P_elec_min");
	TParameter<Double_t>* param_Theta_posi_min = (TParameter<Double_t>*) file_in->Get("Theta_posi_min");
	TParameter<Double_t>* param_Theta_posi_max = (TParameter<Double_t>*) file_in->Get("Theta_posi_max");
	TParameter<Double_t>* param_P_posi_min = (TParameter<Double_t>*) file_in->Get("P_posi_min");
	if (param_seed != nullptr) {
		std::cout << "seed           " << param_seed->GetVal() << std::endl;
	}
	if (param_Nsim != nullptr) {
		std::cout << "Nsim           " << param_Nsim->GetVal() << std::endl;
	}
	if (param_target != nullptr) {
		std::cout << "target         " << param_target->GetVal() << std::endl;
	}
	if (param_Eb != nullptr) {
		std::cout << "Eb             " << param_Eb->GetVal() << std::endl;
	}
	if (param_tMin != nullptr) {
		std::cout << "tMin           " << param_tMin->GetVal() << std::endl;
	}
	if (param_tMax != nullptr) {
		std::cout << "tMax           " << param_tMax->GetVal() << std::endl;
	}
	if (param_EgMin != nullptr) {
		std::cout << "EgMin          " << param_EgMin->GetVal() << std::endl;
	}
	if (param_EgMax != nullptr) {
		std::cout << "EgMax          " << param_EgMax->GetVal() << std::endl;
	}
	if (param_MinvMin != nullptr) {
		std::cout << "MinvMin        " << param_MinvMin->GetVal() << std::endl;
	}
	if (param_MinvMax != nullptr) {
		std::cout << "MinvMax        " << param_MinvMax->GetVal() << std::endl;
	}
	if (param_q2Cut != nullptr) {
		std::cout << "q2Cut          " << param_q2Cut->GetVal() << std::endl;
	}
	if (param_target_l != nullptr) {
		std::cout << "target_l       " << param_target_l->GetVal() << std::endl;
	}
	if (param_target_off != nullptr) {
		std::cout << "target_off     " << param_target_off->GetVal() << std::endl;
	}
	if (param_Dterm != nullptr) {
		std::cout << "Dterm          " << param_Dterm->GetVal() << std::endl;
	}
	if (param_Theta_elec_min != nullptr) {
		std::cout << "Theta_elec_min " << param_Theta_elec_min->GetVal() << std::endl;
	}
	if (param_Theta_elec_max != nullptr) {
		std::cout << "Theta_elec_max " << param_Theta_elec_max->GetVal() << std::endl;
	}
	if (param_P_elec_min != nullptr) {
		std::cout << "P_elec_min     " << param_P_elec_min->GetVal() << std::endl;
	}
	if (param_Theta_posi_min != nullptr) {
		std::cout << "Theta_posi_min " << param_Theta_posi_min->GetVal() << std::endl;
	}
	if (param_Theta_posi_max != nullptr) {
		std::cout << "Theta_posi_max " << param_Theta_posi_max->GetVal() << std::endl;
	}
	if (param_P_posi_min != nullptr) {
		std::cout << "P_posi_min     " << param_P_posi_min->GetVal() << std::endl;
	}

	return 0;
}

