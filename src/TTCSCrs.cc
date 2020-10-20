/* 
 * File:   TTCSCrs.cc
 * Author: rafopar
 * 
 * Created on January 11, 2020, 3:08 PM
 */

#include <TF1.h>
#include <TF2.h>
#include "TTCSCrs.h"
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <KinFunctions.h>

using namespace KinFuncs;

static double sq(double x) {
	return x * x;
}

TTCSCrs::TTCSCrs() {
	TString dat="CFFs_DD_Feb2012.dat";
	if (gSystem->AccessPathName(dat))
		dat=gSystem->Getenv("TCSGEN_DIR")+TString("/")+dat;

	// FIXME
	f_BH = new TF2("f_BH", BH_crs_section, 0, 360, 0, 180, 7);
	f_INT = new TF2("f_INT", INT_crs_section, 0., 360., 0., 180., 12);
	gp = new GPDs(dat, 17, 17, 9, 1.49, -0.20, 0.072);
}

TTCSCrs::TTCSCrs(double a_s, double a_Q2, double a_t) {

	TString dat="CFFs_DD_Feb2012.dat";
	if (gSystem->AccessPathName(dat))
		dat=gSystem->Getenv("TCSGEN_DIR")+TString("/")+dat;

	Set_SQ2t(a_s, a_Q2, a_t);
	iM1 = M_p;
	iM2 = M_p;
	iweight = -1;
	// FIXME
	f_BH = new TF2("f_BH", BH_crs_section, 0, 360, 0, 180, 7);
	f_BH->SetParameters(is, iQ2, it, iM1, iM2);

	f_INT = new TF2("f_INT", INT_crs_section, 0., 360., 0., 180., 12);
	double eta = iQ2 / (2 * (is - M_p * M_p) - iQ2);
	gp = new GPDs(dat, 17, 17, 9, iQ2, it, eta);
	gp->Set_q2_t_eta(iQ2, it, eta);
}

TTCSCrs::~TTCSCrs() {
	delete f_BH;
	delete f_INT;
	delete gp;
}

double TTCSCrs::BH_crs_section( double *x, double *par)
{
	double phi = x[0]/radian;
	double theta = x[1]/radian;

	double s = par[0];
	double Q2 = par[1];
	double t = par[2];
	double M_p1 = par[3];
	double M_p2 = par[4];
	double iweight = par[5];
	// FIXME
	double b = par[6];

	double M2_1 = sq(M_p1);
	double M2_2 = sq(M_p2);
	double M2_avg = 0.5*(M2_1 + M2_2);

	double weight;

	double beta = sqrt( 1 - (4*m_e*m_e)/Q2 );
	double r = sqrt( sq(s - Q2 - M2_2) - 4*Q2*M2_2 );
	// From Byckling equation (4.4.9).
	double cos_TH_Cm = (2*s*(t - M2_1 - M2_2) + (s + M2_1)*(s + M2_2 - Q2)) / sqrt(Lambda(s,M2_1,0)*Lambda(s,M2_2,Q2));
	double sin_TH_Cm = sqrt( 1 - cos_TH_Cm*cos_TH_Cm );
	double Delta_Perp = sin_TH_Cm*r/(2*sqrt(s));
	double a = beta*r*cos(theta);

	// TODO: Fix this too.
	/*
	double b = beta*((Q2*(s - M2_1 - Q2) + t*(s - M2_1 + Q2))/r)*cos(theta) -
		beta*( (2*(s - M2_1)*sqrt(Q2)*Delta_Perp)/r )*sin(theta)*cos(phi);
	*/

	double L_BH = (sq(Q2 - t) - b*b)/4.;
	double L0_BH = Q2*Q2*sin(theta)*sin(theta)/4.;
	// TODO: Fix this too.
	double A_BH = sq((s - M2_1)*Delta_Perp) - t*a*(a + b) - M2_avg*b*b - t*(4*M2_avg - t)*Q2 +
		(m_e*m_e/L_BH)*(sq((Q2 - t)*(a + b) - (s - M2_1)*b) + t*(4*M2_avg - t)*sq(Q2 - t));

	double B_BH = sq(Q2 + t) + b*b + 8*m_e*m_e*Q2 - 4*m_e*m_e*(t + 2*m_e*m_e)*sq(Q2 - t)/L_BH;

	// TODO: Fix this too.
	double F1p = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M2_avg) ))*( 1 - 2.79*t/( 4*M2_avg ) );
	double F2p = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M2_avg) ))*(ammp - 1);

	if( iweight == 1 )
	{ weight = L_BH/L0_BH; }
	else
	{
		weight = 1.;
	}

	// TODO: Fix this too.
	double crs = (1./(2*PI))*weight*sin(theta)*(TMath::Power(alpha_em, 3)/(4*PI*sq(s - M2_1)) )*(beta/(-t*L_BH))*
		( (A_BH/(-t))*(sq(F1p) - (t/(4*M2_avg))*sq(F2p) ) + sq(F1p + F2p)*B_BH/2 )*
		0.389379*1e9;

	return crs;
}

//==================The Interference term ======================

double TTCSCrs::INT_crs_section(double *x, double *par) {
	double phi = x[0] / radian;
	double theta = x[1] / radian;

	double s = par[0];
	double Q2 = par[1];
	double t = par[2];
	double sc_D = par[3]; // Scale of Dterm;
	double weight;

	double ImH = par[5];
	double ReH = par[6];
	double ImE = par[7];
	double ReE = par[8];
	double ImHtild = par[9];
	double ReHtild = par[10];
	double Dterm = par[11];

	double sigma = 1.;

	double beta = sqrt(1 - (4 * m_e * m_e) / Q2);
	double r = sqrt((s - Q2 - M_p * M_p)*(s - Q2 - M_p * M_p) - 4 * Q2 * M_p * M_p);
	double tau = Q2 / (s - M_p * M_p);
	double eta = tau / (2 - tau);
	double cos_TH_Cm = (2 * s * (t - 2 * M_p * M_p) + (s + M_p * M_p)*(s + M_p * M_p - Q2)) / sqrt(Lambda(s, M_p*M_p, 0) * Lambda(s, M_p*M_p, Q2));
	double sin_TH_Cm = sqrt(1 - cos_TH_Cm * cos_TH_Cm);
	double Delta_Perp = sin_TH_Cm * r / (2 * sqrt(s));
	//double Delta_Perp = sqrt((-t)*(1 - tau) - tau*tau*M_p*M_p );
	double a = beta * r * cos(theta);
	double b = sigma * beta * sqrt((Q2 - t)*(Q2 - t) - TMath::Power((2 * (s - M_p * M_p) * sqrt(Q2) * Delta_Perp) / r, 2)) * cos(theta) -
		beta * ((2 * (s - M_p * M_p) * sqrt(Q2) * Delta_Perp) / r) * sin(theta) * cos(phi);
	double L_BH = ((Q2 - t)*(Q2 - t) - b * b) / 4.;
	double L0_BH = Q2 * Q2 * sin(theta) * sin(theta) / 4.;

	double F1p = (1. / ((1 - t / 0.71)*(1 - t / 0.71)))*(1 / (1 - t / (4. * M_p * M_p)))*(1 - 2.79 * t / (4 * M_p * M_p));
	double F2p = (1. / ((1 - t / 0.71)*(1 - t / 0.71)))*(1 / (1 - t / (4. * M_p * M_p)))*(ammp - 1);
	//double F1n = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*(-ammn*t/(4*M_p*M_p));
	//double F2n = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*ammn;

	double t_min = -4 * eta * eta * M_p * M_p / (1 - eta * eta);

	double M2_int = 2 * sqrt(t_min - t) / M_p * (1 - eta) / (1 + eta)*(F1p * (ReH + sc_D * Dterm) - eta * (F1p + F2p) * ReHtild -
			(t / (4 * M_p * M_p)) * F2p * (ReE - sc_D * Dterm));

	if (par[4] == 1) {
		weight = L_BH / L0_BH;
	} else {
		weight = 1;
	}
	// TODO: Check 1/2pi factor.
	double sigma_int = (1/(2*PI))*sin(theta) * weight * (-TMath::Power(alpha_em, 3)) / (4 * PI * s * s)*(-1 / t)*(M_p / sqrt(Q2))*(1 / (tau * sqrt(1 - tau))) * L0_BH / L_BH *
		//double sigma_int = weight*(-TMath::Power(alpha_em, 3))/(4*PI*s*s)*(-1/t)*(M_p/sqrt(Q2))*(1/(tau*sqrt(1 - tau)))*L0_BH/L_BH*
		cos(phi)*(1 + cos(theta) * cos(theta)) / sin(theta) * M2_int *
		0.389379 * 1e9;

	return sigma_int; // This is dsigma/dQ2dts\thetad\phi pbn*GeV-4*rad-1


}

double TTCSCrs::Eval_BH(double a_phi, double a_th) const {
	f_BH->SetParameters(is, iQ2, it, iM1, iM2, iweight, 0);

	return f_BH->Eval(a_phi, a_th) / sin(a_th / radian); // 1/sin(theta) is for getting d\sigma/dcos(\theta)
}

double TTCSCrs::Eval_BH(double a_s, double a_Q2, double a_t, double a_M1, double a_M2, double a_weight, double b, double a_phi, double a_th) const {
	// FIXME
	f_BH->SetParameters(a_s, a_Q2, a_t, a_M1, a_M2, a_weight, b);
	return f_BH->Eval(a_phi, a_th) / sin(a_th / radian); // 1/sin(theta) is for getting d\sigma/dcos(\theta)

}

double TTCSCrs::Eval_INT(double a_phi, double a_th, double a_sc_D) const {

	double eta = iQ2 / (2 * (is - M_p * M_p) - iQ2);
	gp->Set_q2_t_eta(iQ2, it, eta);

	double ImH = gp->GetImH();
	double ReH = gp->GetReH();
	double ImE = gp->GetImE();
	double ReE = gp->GetReE();
	double ImHtild = gp->GetImHtild();
	double ReHtild = gp->GetReHtild();
	double Dterm = gp->GetDterm();

	f_INT->SetParameters(is, iQ2, it, a_sc_D, iweight, ImH, ReH, ImE, ReE, ImHtild, ReHtild);
	f_INT->SetParameter(11, Dterm);
	return f_INT->Eval(a_phi, a_th) / sin(a_th / radian); // 1/sin(theta) is for getting d\sigma/dcos(\theta)
}

double TTCSCrs::Eval_INT(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th, double a_sc_D) const {
	double eta = a_Q2 / (2 * (a_s - M_p * M_p) - a_Q2);
	gp->Set_q2_t_eta(a_Q2, a_t, eta);

	double ImH = gp->GetImH();
	double ReH = gp->GetReH();
	double ImE = gp->GetImE();
	double ReE = gp->GetReE();
	double ImHtild = gp->GetImHtild();
	double ReHtild = gp->GetReHtild();
	double Dterm = gp->GetDterm();
	//  cout<<"ImH   ReH  ImE  ReE  ImHtild ReHtild Dterm  "<<ImH<<"   "<<ReH<<"   "<<ImE<<"   "<<ReE<<"   "<<ImHtild<<"   "<<ReHtild<<"   "<<Dterm<<endl;

	f_INT->SetParameters(a_s, a_Q2, a_t, a_sc_D, a_weight, ImH, ReH, ImE, ReE, ImHtild, ReHtild);
	f_INT->SetParameter(11, Dterm);
	return f_INT->Eval(a_phi, a_th) / sin(a_th / radian); // 1/sin(theta) is for getting d\sigma/dcos(\theta)
}

void TTCSCrs::Set_SQ2t(double a_s, double a_Q2, double a_t) {
	is = a_s;
	it = a_t;
	iQ2 = a_Q2;
}

void TTCSCrs::Set_M1M2(double a_M1, double a_M2) {
	iM1 = a_M1;
	iM2 = a_M2;
}

double TTCSCrs::Integral_BH_phi_th(double a_phi_min, double a_phi_max, double a_th_min, double a_th_max) {

	// FIXME
	f_BH->SetParameters(is, iQ2, it, iM1, iM2, iweight, 0);
	return f_BH->Integral(a_phi_min, a_phi_max, a_th_min, a_th_max);
}

void TTCSCrs::Set_Weight(double a_weight) {
	iweight = a_weight;
}

void TTCSCrs::Set_sc_D(double a_sc_D) {
	isc_D = a_sc_D;
}

void TTCSCrs::Draw_BH(const char* option) {
	f_BH->SetNpx(900);
	f_BH->SetNpy(900);
	// FIXME
	f_BH->SetParameters(is, iQ2, it, iM1, iM2, iweight, 0);
	f_BH->Draw(option);
}

void TTCSCrs::Draw_INT(const char* option, double a_sc_D) {
	double eta = iQ2 / (2 * (is - M_p * M_p) - iQ2);
	gp->Set_q2_t_eta(iQ2, it, eta);

	double ImH = gp->GetImH();
	double ReH = gp->GetReH();
	double ImE = gp->GetImE();
	double ReE = gp->GetReE();
	double ImHtild = gp->GetImHtild();
	double ReHtild = gp->GetReHtild();
	double Dterm = gp->GetDterm();

	f_INT->SetParameters(is, iQ2, it, a_sc_D, iweight, ImH, ReH, ImE, ReE, ImHtild, ReHtild);
	f_INT->SetParameter(11, Dterm);
	f_INT->SetNpx(500);
	f_INT->SetNpy(500);
	f_INT->Draw(option);
}

TH2D *TTCSCrs::Get_BH_Crs_Histogream_ThPhi(const char *name, int weight) {
	TCanvas *ctmp = new TCanvas();
	// FIXME
	f_BH->SetParameters(is, iQ2, it, iM1, iM2, weight, 0);
	f_BH->SetNpx(500);
	f_BH->SetNpy(500);
	f_BH->Draw("colz");
	TH2D *h_tmp = (TH2D*) f_BH->GetHistogram();
	h_tmp->SetName(name);
	delete ctmp;
	return h_tmp;
}

TH2D *TTCSCrs::Get_INT_Crs_Histogream_ThPhi(const char *name, int weight) {
	TCanvas *ctmp = new TCanvas(); // Without this I was getting crazy results in utility programs :-)
	double eta = iQ2 / (2 * (is - M_p * M_p) - iQ2);
	gp->Set_q2_t_eta(iQ2, it, eta);

	double ImH = gp->GetImH();
	double ReH = gp->GetReH();
	double ImE = gp->GetImE();
	double ReE = gp->GetReE();
	double ImHtild = gp->GetImHtild();
	double ReHtild = gp->GetReHtild();
	double Dterm = gp->GetDterm();

	f_INT->SetParameters(is, iQ2, it, isc_D, weight, ImH, ReH, ImE, ReE, ImHtild, ReHtild);
	f_INT->SetParameter(11, Dterm);
	f_INT->SetNpx(500);
	f_INT->SetNpy(500);
	f_INT->Draw("colz");
	TH2D *h_tmp = (TH2D*) f_INT->GetHistogram();
	h_tmp->SetName(name);
	delete ctmp;
	return h_tmp;
}

