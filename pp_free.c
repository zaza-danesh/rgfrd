//
//  pp_free.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/06/13.
//
//

 
#include "pp_free.h"

double pp_free(double theta, double r, double t, phys_parameters_bimol1 physparams){

	double result = 0.0;
	const double ka = physparams.ka;
	const double D = physparams.D;
	const double sigma = physparams.sigma;
	const double r0 = physparams.r0;

	const double Dt = D*t;
	const double Dt2 = 2*Dt;
	const double Dt4 = 2*Dt2;
	const double Dt2_r = 1.0/Dt2;
	const double Dt4_r = Dt2_r/2.0;
	const double r_sqr = r*r;
	const double r0_sqr = r0*r0;
	const double r_r0 = r*r0;
	const double cos_theta = cos(theta);
	const double Pi_cube = Pi*Pi*Pi;
	const double r_r0_over_Dt2 = r_r0*Dt2_r;
	const double r_sqr_p_r0_sqr = r_sqr + r0_sqr;
	const double r_sqr_p_r0_sqr_over_Dt4 = r_sqr_p_r0_sqr*Dt4_r;

	

	const double term1 = expm1(r_r0_over_Dt2 - r_sqr_p_r0_sqr_over_Dt4);
	const double term2 = expm1(r_r0_over_Dt2*cos_theta - r_sqr_p_r0_sqr_over_Dt4);

	const double num = term1 - term2;
	const double den = 4*r_r0*sqrt(Pi_cube*Dt);

	result = num/den;

	return result;
}







