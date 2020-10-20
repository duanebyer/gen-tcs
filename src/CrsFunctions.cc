#include <CrsFunctions.h>
#include <KinFunctions.h>
#include <iostream>
#include <cmath>

using namespace KinFuncs;

static double sq(double x) {
	return x * x;
}

double CrsFuncs::BH(
		double s, double t,
		double p2, double p2prime,
		double q2prime, double k2,
		double a, double b,
		int weight) {
	// Fine structure constant.
	double alpha_em = 1. / 137;
	// Magnetic moment of the proton.
    double ammp = 2.793;
	double p2avg = 0.5*(p2 + p2prime);
	double beta = sqrt(1 - (4*k2)/q2prime);
	double r = sqrt(sq(s - q2prime - p2prime) - 4*q2prime*p2prime);
	// From Byckling equation (4.4.9).
	double cos_TH_Cm = (2*s*(t - p2 - p2prime) + (s + p2)*(s + p2prime - q2prime))
		/ sqrt(Lambda(s,p2,0)*Lambda(s,p2prime,q2prime));
	double sin_TH_Cm = sqrt(1 - cos_TH_Cm*cos_TH_Cm);
	double Delta_Perp = sin_TH_Cm*r/(2*sqrt(s));

	double L_BH = (sq(q2prime - t) - sq(b))/4;
	// TODO: Simplify sin(acos(...)) expression in L0_BH.
	double theta = acos(a / (beta * r));
	double L0_BH = sq(q2prime*sin(theta))/4;
	// TODO: Fix this too.
	double A_BH =
		sq((s - p2)*Delta_Perp)
		- t*a*(a + b)
		- 0.5*(p2 + p2prime)*b*b
		- t*(2*p2 + 2*p2prime - t)*q2prime
		+ (k2/L_BH)*(
			sq((q2prime - t)*(a + b) - (s - p2)*b)
			+ t*(2*p2 + 2*p2prime - t)*sq(q2prime - t));

	double B_BH =
		sq(q2prime + t) + b*b + 8*k2*q2prime - 4*k2*(t + 2*k2)*sq(q2prime - t)/L_BH;

	// TODO: Test if cross section changes if proton mass or p2avg is used for
	// form factor.
	double F1p = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*p2avg) ))*( 1 - 2.79*t/( 4*p2avg) );
	double F2p = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*p2avg) ))*(ammp - 1);

	double weight_factor = weight == 1 ? L_BH / L0_BH : 1;

	double crs =
		(1/(2*M_PI))*weight_factor
		*(alpha_em*alpha_em*alpha_em/(4*M_PI*sq(s - p2)))
		*(beta/(-t*L_BH))
		*((A_BH/(-t))*(sq(F1p) - (t/(4*p2avg))*sq(F2p)) + sq(F1p + F2p)*B_BH/2)
		*0.389379*1e9;

	return crs;
}

double CrsFuncs::BHInt(
		double s, double t,
		double p2, double p2prime,
		double q2prime, double k2,
		double a, double b,
		double Dterm,
		int weight) {
	return 0;
}

