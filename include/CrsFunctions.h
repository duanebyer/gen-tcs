#ifndef _CRS_FUNCTIONS_H_
#define _CRS_FUNCTIONS_H_

namespace CrsFuncs {

	double BH(
		// Kinematic variables.
		double s, double t,
		double p2, double p2prime,
		double q2prime, double k2,
		// Kinematic variables related to angle.
		double a, double b,
		// Kind of weighting to use (zero or one).
		int weight);

	double BHInt(
		double s, double t,
		double p2, double p2prime,
		double q2prime, double k2,
		double a, double b,
		double Dterm,
		int weight);
}

#endif

