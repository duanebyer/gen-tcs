/* 
 * File:   TCSGen.cc
 * Author: rafopar
 *
 * Created on January 11, 2020, 2:53 PM
 */

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <KinFunctions.h>
#include <CrsFunctions.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParameter.h>
#include <TRandom2.h>
#include <TTCSCrs.h>
#include <TTCSKine.h>
#include <TTree.h>

const double PI = TMath::Pi();
const double PROTON_MASS = 0.9383;
const double NEUTRON_MASS = 0.9396;
const double ELECTRON_MASS = 0.00051;
const double DEUTERON_MASS = 1.8756;

enum class Target {
	PROTON,
	DEUTERON
};

double Sq(double x) {
	return x * x;
}

int main(int argc, char** argv) {
	std::ifstream config_file("GenOptions.dat");
	std::map<std::string, std::string> m_Settings;
	if (config_file.is_open()) {
		while (!config_file.eof()) {
			std::string Key;
			std::string Val;
			config_file >> Key;
			config_file >> Val;
			m_Settings[Key] = Val;
		}
	} else {
		std::cout << "Failed to open GenOptions.dat." << std::endl;
		return 1;
	}

	int seed              = std::stoi(m_Settings.at("seed"));
	int Nsim              = std::stoi(m_Settings.at("Nsim"));
	double Eb_user        = std::stod(m_Settings.at("Eb"));

	// Limits given by the user. These will be used if possible, but limits from
	// kinematics will take precendence.
	double t_minUser      = std::stod(m_Settings.at("tMin"));
	double t_maxUser      = std::stod(m_Settings.at("tMax"));
	double q2_cut         = std::stod(m_Settings.at("q2Cut"));
	double Minv_MinUser   = std::stod(m_Settings.at("MinvMin"));
	double Minv_MaxUser   = std::stod(m_Settings.at("MinvMax"));
	// TODO: Eg is measured in the frame where the target is at rest. Determine
	// if this is appropriate.
	double Eg_minUser     = std::stod(m_Settings.at("EgMin"));
	double Eg_maxUser     = std::stod(m_Settings.at("EgMax"));

	// Treat target as deuteron.
	Target target_type = Target::PROTON;
	if (m_Settings.at("target") == "proton" || m_Settings.at("target") == "0") {
		target_type = Target::PROTON;
	} else if (m_Settings.at("target") == "deuteron" || m_Settings.at("target") == "1") {
		target_type = Target::DEUTERON;
	} else {
		std::cout << "Invalid target type \'" << m_Settings.at("target") << "\'." << std::endl;
		return 1;
	}

	// Target parameters.
	double l_targ         = std::stod(m_Settings.at("target_l"));
	double offset_targ    = std::stod(m_Settings.at("target_off"));

	// D-term.
	double Dterm          = std::stod(m_Settings.at("Dterm"));

	// Limits on the particle angles.
	double Theta_elec_min = std::stod(m_Settings.at("Theta_elec_min"));
	double Theta_elec_max = std::stod(m_Settings.at("Theta_elec_max"));
	double P_elec_min     = std::stod(m_Settings.at("P_elec_min"));
	double Theta_posi_min = std::stod(m_Settings.at("Theta_posi_min"));;
	double Theta_posi_max = std::stod(m_Settings.at("Theta_posi_max"));;
	double P_posi_min     = std::stod(m_Settings.at("P_posi_min"));;

	// Parameters to be stored in the file for later reference.
	TParameter<Int_t> param_seed("seed", seed);
	TParameter<Int_t> param_Nsim("Nsim", Nsim);
	TParameter<Int_t> param_target("target", (int) target_type);
	TParameter<Double_t> param_Eb("Eb", Eb_user);
	TParameter<Double_t> param_tMin("tMin", t_minUser);
	TParameter<Double_t> param_tMax("tMax", t_maxUser);
	TParameter<Double_t> param_EgMin("EgMin", Eg_minUser);
	TParameter<Double_t> param_EgMax("EgMax", Eg_maxUser);
	TParameter<Double_t> param_MinvMax("MinvMax", Minv_MaxUser);
	TParameter<Double_t> param_MinvMin("MinvMin", Minv_MinUser);
	TParameter<Double_t> param_q2Cut("q2Cut", q2_cut);
	TParameter<Double_t> param_target_l("target_l", l_targ);
	TParameter<Double_t> param_target_off("target_off", offset_targ);
	TParameter<Double_t> param_Dterm("Dterm", Dterm);
	TParameter<Double_t> param_Theta_elec_min("Theta_elec_min", Theta_elec_min);
	TParameter<Double_t> param_Theta_elec_max("Theta_elec_max", Theta_elec_max);
	TParameter<Double_t> param_P_elec_min("P_elec_min", P_elec_min);
	TParameter<Double_t> param_Theta_posi_min("Theta_posi_min", Theta_posi_min);
	TParameter<Double_t> param_Theta_posi_max("Theta_posi_max", Theta_posi_max);
	TParameter<Double_t> param_P_posi_min("P_posi_min", P_posi_min);

	TRandom2 rand;
	rand.SetSeed(seed);
	std::cout << "Random seed: " << seed << std::endl;
	std::cout << "Test random number: " << rand.Uniform(0, 1) << std::endl;

	TFile* file_out = new TFile("tcs-gen.root", "Recreate");

	TH2D *h_ph_h_ph_cm1 = new TH2D("h_ph_h_ph_cm1", "", 200, 0., 360., 200, 0., 360.);
	TH2D *h_th_g_th_cm1 = new TH2D("h_th_g_th_cm1", "", 200, 0., 180., 200, 0., 180.);

	//================= Definition of Tree Variables =================
	double Eg;
	double t;
	double Q2;
	double psf, crs_BH, crs_INT;
	double flux_factor;
	double flux_brem;
	TLorentzVector L_em, L_ep, L_p, L_pprime, L_g, L_gprime, L_k, L_kprime;

	TTree* events = new TTree("tr1", "TCS MC events");
	events->Branch("L_em", "TLorentzVector", &L_em, 3200);
	events->Branch("L_ep", "TLorentzVector", &L_ep, 3200);
	events->Branch("L_g", "TLorentzVector", &L_g, 3200);
	events->Branch("L_k", "TLorentzVector", &L_kprime, 3200);
	events->Branch("L_k_i", "TLorentzVector", &L_k, 3200);
	events->Branch("L_prot", "TLorentzVector", &L_pprime, 3200);
	events->Branch("L_prot_i", "TLorentzVector", &L_p, 3200);
	events->Branch("Eg", &Eg, "Eg/D");
	events->Branch("Q2", &Q2, "Q2/D");
	events->Branch("t", &t, "t/D");
	events->Branch("psf", &psf, "psf/D");
	events->Branch("flux_factor", &flux_factor, "flux_factor/D");
	events->Branch("flux_brem", &flux_brem, "flux_brem/D");
	events->Branch("crs_BH", &crs_BH, "crs_BH/D");
	events->Branch("crs_INT", &crs_INT, "crs_INT/D");

	int event_idx=0;

	double a = 0.0456; // GeV
	double b = 0.2719; // GeV
	// eq 77
	TF1* deuteron_k = new TF1("deuteron_k", "x^2*(1/(x^2+[0]^2)-1/(x^2+[1]^2))^2");
	deuteron_k->SetNpx(400);
	deuteron_k->SetParameters(0, a);
	deuteron_k->SetParameters(1, b);

	TTCSCrs crs_lmlp;
	int Nsim_tot = 0;
	while (event_idx < Nsim) {
		Nsim_tot += 1;

		if (target_type == Target::PROTON) {
			L_p.SetPxPyPzE(0., 0., 0., PROTON_MASS);
		} else if (target_type == Target::DEUTERON) {
			double k_max = 1. / (2.*DEUTERON_MASS)*(Sq(DEUTERON_MASS) - Sq(NEUTRON_MASS)); // GeV
			//double k_max = 0.2;
			double k_mag = deuteron_k->GetRandom(0, k_max);
			double kx, ky, kz;
			rand.Sphere(kx, ky, kz, k_mag);
			// The proton is off-shell, but the neutron created when the
			// deuteron splits up is on-shell, so we can use that to figure out
			// the proton energy.
			double energy = DEUTERON_MASS - TMath::Sqrt(Sq(k_mag) + Sq(NEUTRON_MASS));
			L_p.SetPxPyPzE(-kx, -ky, -kz, energy);
		} else {
			std::cout << "Cannot use target type " << (int) target_type << std::endl;
			return 1;
		}

		// Set electron beam.
		L_k.SetPxPyPzE(0., 0., TMath::Sqrt(Sq(Eb_user) - Sq(ELECTRON_MASS)), Eb_user);

		// Transform into target frame.
		TVector3 target_boost_inv = L_p.BoostVector();
		TVector3 target_boost = -target_boost_inv;
		TRotation target_rotation, target_rotation_inv;
		L_k.Boost(target_boost);
		L_p.Boost(target_boost);
		target_rotation_inv.SetZAxis(L_k.Vect());
		target_rotation = target_rotation_inv.Inverse();
		L_k.Transform(target_rotation);
		L_p.Transform(target_rotation);
		double Eb = L_k.E();
		double Mp = L_p.M();
		double Mpprime = PROTON_MASS;
		double Mp2 = Mp*Mp;
		double Mpprime2 = Mpprime*Mpprime;

		// Find the direction of the photon in the lab frame.
		TLorentzVector L_g_dir(0., 0., 1., 1.);
		L_g_dir.Transform(target_rotation_inv);
		L_g_dir.Boost(target_boost_inv);
		// Then set the energy bounds in the lab frame and transform back.
		TLorentzVector L_g_minUser = Eg_minUser * L_g_dir;
		TLorentzVector L_g_maxUser = Eg_maxUser * L_g_dir;
		L_g_minUser.Boost(target_boost);
		L_g_maxUser.Boost(target_boost);
		L_g_minUser.Transform(target_rotation);
		L_g_maxUser.Transform(target_rotation);
		// Find the new energy bounds in the target frame.
		Eg_minUser = L_g_minUser.E();
		Eg_maxUser = L_g_maxUser.E();

		// The beam energy provides the kinematic limit on Eg.
		double Eg_minKine = 0.;
		double Eg_min = TMath::Max(Eg_minKine, Eg_minUser);
		double Eg_maxKine = Eb;
		double Eg_max = TMath::Min(Eg_maxKine, Eg_maxUser);
		double psf_Eg = Eg_max - Eg_min;
		Eg = rand.Uniform(Eg_min, Eg_max);

		// COM energy: s = (p + q)^2
		double s = Mp2 + 2*Mp*Eg;

		// Check that there is enough energy in the COM frame to create the
		// leptons and the final proton.
		if (s < Sq(Minv_MinUser + Mpprime)) {
			continue;
		}

		// Check that the photon energy is above the threshold.
		// TODO: This is the same as the previous condition on s, so why is it
		// here a second time?
		double Eg_threshold = (Sq(Minv_MinUser + Mpprime) - Mp2) / (2*Mp);
		if (Eg < Eg_threshold) {
			continue;
		}

		L_g.SetPxPyPzE(0, 0, Eg, Eg);
		L_kprime = L_k - L_g;

		// We will mostly use Q2 instead of Minv.
		double Q2MinUser = Minv_MinUser*Minv_MinUser;

		// TODO: Check whether Q2MinUser is still the appropriate choice if it
		// is possible for Mp != Mpprime. Alternative would be to take the min/
		// max of Q2MinUser and Q2MaxUser.
		// Addenum: I think this is valid, because Q2 is chosen after this, and
		// the upper limit on Q2 is determined by t. So, Q2min determines
		// allowed t range, which in turn determines Q2max. Need to read
		// Byckling to be sure.
		double t_min = TMath::Max(
			KinFuncs::T_max(0., Mp2, Q2MinUser, Mpprime2, s), t_minUser);
		double t_max = TMath::Min(
			KinFuncs::T_min(0., Mp2, Q2MinUser, Mpprime2, s), t_maxUser);
		double psf_t = t_max - t_min;
		if (t_min > t_max) {
			continue;
		}
		t = rand.Uniform(t_min, t_max);

		// For a given Q2MinUser, we should check the EgMinUser. If it is below
		// the threshold to produce Minv_MinUser, then Eg_min should be adjusted
		// to the threshold for Minv production.

		// TODO: Check that tcs_kin1 works correctly with off-shell proton mass.
		TTCSKine tcs_kin1(Mp, Eb);

		flux_factor = KinFuncs::N_EPA(Eb, Eg, q2_cut);
		flux_brem = KinFuncs::Brem_Approx(Eg, Eb, l_targ);
		// The flux factor can sometimes be below zero due to numerical issues.
		if (flux_factor < 0.0) {
			flux_factor = 0.0;
		}

		// Equation for s2 p121(5/5.11) in Byckling
		// The kinematic upper bound on Q2.
		double Q2maxKine =
			s + Mpprime2 - (1/(2*Mp2))*((s + Mp2)*(Mp2 + Mpprime2 - t)
			- TMath::Sqrt(KinFuncs::Lambda(s, Mp2, 0.)*KinFuncs::Lambda(t, Mp2, Mpprime2)));
		double Q2min = Q2MinUser;
		double Q2max = TMath::Min(Q2maxKine, Sq(Minv_MaxUser));
		double psf_Q2 = Q2max - Q2min;
		if (Q2min > Q2max) {
			continue;
		}
		Q2 = rand.Uniform(Q2min, Q2max);

		// Byukling Kayanti (4.9)
		// Polar angle in COM frame.
		double th_qprime = acos((2*s*(t - Mp2 - Mpprime2) + (s + Mpprime2 - Q2)*(s + Mp2)) / TMath::Sqrt(KinFuncs::Lambda(s, 0, Mp2) * KinFuncs::Lambda(s, Q2, Mpprime2)));
		double th_pprime = PI + th_qprime;
		// COM momentum (same for p' and q').
		double Pprime = 0.5 * TMath::Sqrt(KinFuncs::Lambda(s, Q2, Mpprime2) / s);

		// cos(th) is in range (-1, 1).
		double psf_cos_th = 2.; 
		double psf_phi_cm = 2 * PI;

		double cos_th = rand.Uniform(-1., 1.);
		double sin_th = TMath::Sqrt(1 - cos_th * cos_th);
		double phi_cm = rand.Uniform(0., 2. * PI);

		// Energy of lepton in the rest frame of qprime.
		double El = 0.5 * sqrt(Q2);
		double Pl = sqrt(Sq(El) - Sq(ELECTRON_MASS));

		// Move to COM frame from target frame.
		TLorentzVector L_cm;
		L_cm.SetPxPyPzE(0., 0., Eg, Mp + Eg);

		// Find final proton and photon momentum in COM frame.
		L_pprime.SetPxPyPzE(
			Pprime * TMath::Sin(th_pprime),
			0.,
			Pprime * TMath::Cos(th_pprime),
			TMath::Sqrt(Pprime * Pprime + Mpprime2));
		L_gprime.SetPxPyPzE(
			Pprime * sin(th_qprime),
			0.,
			Pprime * TMath::Cos(th_qprime),
			TMath::Sqrt(Pprime * Pprime + Q2));

		// Create the leptons in the lepton frame.
		L_em.SetPxPyPzE(
			Pl * sin_th * TMath::Cos(phi_cm),
			Pl * sin_th * TMath::Sin(phi_cm),
			Pl * cos_th,
			El);
		L_ep.SetPxPyPzE(
			-Pl * sin_th * TMath::Cos(phi_cm),
			-Pl * sin_th * TMath::Sin(phi_cm),
			-Pl*cos_th,
			El);

		// Rotate in order to get Z axis be antiparallel to the p_prime
		// direction in the CM frame
		L_em.RotateY(th_qprime); 
		L_ep.RotateY(th_qprime);

		// Move to the CM frame from the lepton frame.
		L_em.Boost(L_gprime.BoostVector());
		L_ep.Boost(L_gprime.BoostVector());

		// Move to the target frame from the CM frame.
		L_em.Boost(L_cm.BoostVector());
		L_ep.Boost(L_cm.BoostVector());
		L_gprime.Boost(L_cm.BoostVector());
		L_pprime.Boost(L_cm.BoostVector());

		double psf_phi_lab = 2 * PI;
		double phi_rot = rand.Uniform(0., 2. * PI);

		L_pprime.RotateZ(phi_rot);
		L_gprime.RotateZ(phi_rot);
		L_em.RotateZ(phi_rot);
		L_ep.RotateZ(phi_rot);
		tcs_kin1.SetLemLepLp(L_em, L_ep, L_pprime);

		h_ph_h_ph_cm1->Fill(phi_cm * TMath::RadToDeg(), tcs_kin1.GetPhi_cm());
		h_th_g_th_cm1->Fill(TMath::ACos(cos_th) * TMath::RadToDeg(), tcs_kin1.GetTheta_cm());

		psf = psf_Eg * psf_t * psf_Q2 * psf_phi_lab * psf_cos_th * psf_phi_cm;

		double a = 2*(L_pprime).Dot(L_em - L_ep);
		double b = 2*(L_p - L_pprime).Dot(L_em - L_ep);
		//crs_BH = crs_lmlp.Eval_BH(s, Q2, t, Mp, Mpprime, -1, b, tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm());
		crs_BH = CrsFuncs::BH(s, t, Mp2, Mpprime2, Q2, Sq(ELECTRON_MASS), a, b, -1);

		// We want to make sure the kinematics is inside the grid of CFFs,
		// otherwise the cross section is not defined.
		double eta = Q2 / (2 * (s - Mp * Mp) - Q2);
		crs_INT = CrsFuncs::BHInt(s, t, Mp, Mpprime, Q2, Sq(ELECTRON_MASS), a, b, Dterm, -1);
		/*if (Q2 < 9. && -t < 0.8 && eta < 0.8 && eta > 0.06) {
			crs_INT = crs_lmlp.Eval_INT(s, Q2, t, -1., tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm(), Dterm);
		} else {
			crs_INT = 0;
		}*/
		
		// Move to the lab frame from target frame.
		// TODO: Should all momenta be included in here?
		L_p.Transform(target_rotation_inv);
		L_pprime.Transform(target_rotation_inv);
		L_em.Transform(target_rotation_inv);
		L_ep.Transform(target_rotation_inv);
		L_k.Transform(target_rotation_inv);
		L_kprime.Transform(target_rotation_inv);
		L_g.Transform(target_rotation_inv);

		L_p.Boost(target_boost_inv);
		L_pprime.Boost(target_boost_inv);
		L_em.Boost(target_boost_inv);
		L_ep.Boost(target_boost_inv);
		L_k.Boost(target_boost_inv);
		L_kprime.Boost(target_boost_inv);
		L_g.Boost(target_boost_inv);

		bool elec_in_bounds =
			L_em.Theta() * TMath::RadToDeg() >= Theta_elec_min
			&& L_em.Theta() * TMath::RadToDeg() < Theta_elec_max
			&& L_em.P() >= P_elec_min;
		bool posi_in_bounds =
			L_ep.Theta() * TMath::RadToDeg() >= Theta_posi_min
			&& L_ep.Theta() * TMath::RadToDeg() < Theta_posi_max
			&& L_ep.P() >= P_posi_min;
		if (!elec_in_bounds || !posi_in_bounds) {
			continue;
		}

		if (!(psf >= 0.0 && crs_BH >= 0.0 && flux_factor >= 0.0 && flux_brem >= 0.0)) {
			std::cout << "CROSS SECTION"              << std::endl;
			std::cout << "  crs_BH:  " << crs_BH      << std::endl;
			std::cout << "FLUX FACTOR"                << std::endl;
			std::cout << "  ff:      " << flux_factor << std::endl;
			std::cout << "  ff_brem: " << flux_brem   << std::endl;
			std::cout << "PSF"                        << std::endl;
			std::cout << "  Eg:      " << psf_Eg      << std::endl;
			std::cout << "  t:       " << psf_t       << std::endl;
			std::cout << "  Q2:      " << psf_Q2      << std::endl;
			std::cout << "  phi_lab: " << psf_phi_lab << std::endl;
			std::cout << "  cos_th:  " << psf_cos_th  << std::endl;
			std::cout << "  phi_cm:  " << psf_phi_cm  << std::endl;
			return 1;
		}

		if ((event_idx % (Nsim / 100) == 0)) {
			std::cout << "Processed " << double(100. * event_idx / double(Nsim)) << "% of events" << std::endl;
		}
		event_idx += 1;

		events->Fill();
	}

	events->Write();
	h_ph_h_ph_cm1->Write();
	h_th_g_th_cm1->Write();

	TParameter<Int_t> param_Nsim_tot("Nsim_tot", Nsim_tot);
	param_seed.Write();
	param_Nsim_tot.Write();
	param_Nsim.Write();
	param_target.Write();
	param_Eb.Write();
	param_tMin.Write();
	param_tMax.Write();
	param_EgMin.Write();
	param_EgMax.Write();
	param_MinvMax.Write();
	param_MinvMin.Write();
	param_q2Cut.Write();
	param_target_l.Write();
	param_target_off.Write();
	param_Dterm.Write();
	param_Theta_elec_min.Write();
	param_Theta_elec_max.Write();
	param_P_elec_min.Write();
	param_Theta_posi_min.Write();
	param_Theta_posi_max.Write();
	param_P_posi_min.Write();

	file_out->Close();

	return 0;
}

