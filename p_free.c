//
//  p_free.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/06/13.
//
//

 
#include "p_free.h"

double p_free(double theta, double r, double t, phys_parameters_bimol1 physparams){

	double result = 0.0;
	const double ka = physparams.ka;
	const double D = physparams.D;
	const double sigma = physparams.sigma;
	const double r0 = physparams.r0;

	const double Dt = D*t;
	const double Dt4 = 4*Dt;
	const double Dt2 = 2*Dt;
	const double r_sqr = r*r;
	const double r0_sqr = r0*r0;
	const double r_r0 = r*r0;
	const double r_r0_r = 1.0/r_r0;
	const double cos_theta = cos(theta);
	const double Pi2 = 2*Pi;
	const double den = sqrt(Pi2*Dt*Dt*Dt)*r_r0_r;

	const double term1 = exp(-(r_sqr + r0_sqr)/Dt4);
	const double term2 = exp (r_r0/Dt2) - exp(r_r0*cos_theta/Dt2);
	double num = term1*term2;

	num = r_sqr*Dt*num;

	result = num/den;
	
	return result;
}







