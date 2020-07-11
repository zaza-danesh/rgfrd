//
//  space_theta_bimol1.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/02/13.
//
 
#include "space_theta_bimol1.h"

double space_theta_bimol1(double theta, double r, double t, double ka, double D, double sigma, double r0){


	// Integration of (p_free + p_corr from 0 to theta to compute the cumulative function)
	struct phys_parameters_bimol1 physparams;
	physparams.ka = ka;
	physparams.D = D;
	physparams.sigma = sigma;
	physparams.r0 = r0;

	const double p_free = pp_free(theta, r, t, physparams);
	const double p_corr = pp_corr(theta,  r,  t,  physparams);;
	const double result =  p_free + p_corr;  


if ( fabs(result) > epsilon){
          return result;
     }
     else{
          return 0.0;
     }


}