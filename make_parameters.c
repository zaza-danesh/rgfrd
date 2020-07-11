//
//  make_parameters.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/25.
//
//

#include "make_parameters.h"

phys_parameters_bimol1 make_parameters_bimol1( double ka, double D, double sigma,  double r0 ){

	struct phys_parameters_bimol1 physparams;

	physparams.ka = ka;
	physparams.D = D;
	physparams.sigma = sigma;
	physparams.r0 = r0;

	return physparams;
}

phys_parameters_bimol2 make_parameters_bimol2(double ka, double D, double sigma, double b, double kd, double gamma_b, double gamma_u){

	struct phys_parameters_bimol2 physparams;

	physparams.ka = ka;
	physparams.D = D;
	physparams.sigma = sigma;
	physparams.b = b;
	physparams.kd = kd;
	physparams.gamma_u = gamma_u;
    physparams.gamma_b = gamma_b;

	return physparams;
}