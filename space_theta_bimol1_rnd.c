//
//  space_theta_bimol1_rnd.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/02/13.
//
//

 
#include "space_theta_bimol1_rnd.h"

double space_theta_bimol1_rnd(double theta, void* params){

struct space_theta_bimol1_parameters *p  = (struct space_theta_bimol1_parameters *) params;

	double result =0.0;
	double ka = p->ka;;
	double D = p->D;
	double sigma = p->sigma;
	double r0 = p->r0;
	double rnd = p->rnd; 
	double t = p->sampled_time;
	double r = p->sampled_r;


	result = space_theta_bimol1(theta, r, t , ka, D, sigma, r0) - rnd;

	return result;
}