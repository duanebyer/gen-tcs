#include <cmath>
#include <KinFunctions.h>

// This is from Byukling Kayanti Formula (6.3)
double KinFuncs::Lambda( double x, double y, double z )
{
  return (x - y - z)*(x - y - z) - 4*y*z;
}

// From Byukling Kayanti Formula (4.5.23)
double KinFuncs::G(double x, double y, double z, double u, double v, double w) {
	return x*y*(x+y) + z*u*(z+u) + v*w*(v+w) + z*w*(x+y) + x*u*v + y*u*w - x*y*(z + u + v + w) - z*u*(x + y + v + w) - v*w*(x + y + z + u);
}

//From Byukling Kayanti Formula (5.14) Page 86
double KinFuncs::T_min( double ma_2, double mb_2, double m1_2, double m2_2, double s) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
{
  return ma_2 + m1_2 - (1/(2*s))*( (s + ma_2 - mb_2)*(s + m1_2 - m2_2) - sqrt( Lambda(s, ma_2, mb_2)*Lambda(s, m1_2, m2_2) ) );
}

//From Byukling Kayanti Formula (5.14) page 86
double KinFuncs::T_max( double ma_2, double mb_2, double m1_2, double m2_2, double s)
{
  return ma_2 + m1_2 - (1/(2*s))*( (s + ma_2 - mb_2)*(s + m1_2 - m2_2) + sqrt( Lambda(s, ma_2, mb_2)*Lambda(s, m1_2, m2_2) ) );
}

double KinFuncs::Q2_min( double s, double Eb, double M )
{
  // M is the target mass;
  double me = 0.00051;
  double Eg = (s - M*M)/(2*M);
  double E_pr = Eb - Eg;
  double P0 = sqrt(Eb*Eb - me*me);
  double P_pr = sqrt(E_pr*E_pr - me*me);
  double Q2min = 2*(Eb*E_pr - P0*P_pr - me*me);

  return Q2min;
}

double KinFuncs::N_EPA(double Eb, double Eg, double Q2_max)
{
  const double alpha = 1./137.;
  const double PI = 3.14159265358979312;

  double x = Eg/Eb;
  double me = 0.00051;
  double Mp = 0.9383;
  double Q2_min = me*me*x*x/(1 - x);
  //return (1/Eb)*alpha/(PI*x)*( (1 - x + x*x/2)*log(Q2_max/Q2_min) - (1 - x));
  return (1./Eg)*alpha/PI*((1. - x + x*x)*log(Q2_max/Q2_min) - (1. - 0.5*x)*(1. - 0.5*x)*log((Eg*Eg + Q2_max)/(Eg*Eg + Q2_min)) - me*me*x*x*(1./Q2_min - 1./Q2_max));
}

double KinFuncs::Brem_Approx(double Eg, double Eb,double l_targ)
{
	double X_0 = 890.;   // LH2 rad. length in cm
	double x = Eg/Eb;
	// 0.5 is for effective luminosity (after integration int_0^l (l-x)dx )
	//return 0.5*(l_targ/X_0)*(4./3.)*(1/Eg-1/Eb);
	return 0.5*(1./Eg)*(l_targ/X_0)*(4./3.)*(1 - x + x*x);
}

