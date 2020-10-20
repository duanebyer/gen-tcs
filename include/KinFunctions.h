/* 
 * File:   KinFunctions.h
 * Author: rafopar
 *
 * Created on January 11, 2020, 4:38 PM
 */

#ifndef _KIN_FUNCTIONS_H
#define _KIN_FUNCTIONS_H

namespace KinFuncs {

    double Lambda(double x, double y, double z);
	double G(double x, double y, double z, double u, double v, double w);
    double T_min(double, double, double, double, double);
    double T_max(double, double, double, double, double);
    double Q2_min(double s, double Eb, double M);
    double N_EPA(double, double, double );
	double Brem_Approx(double, double, double);
}

#endif

