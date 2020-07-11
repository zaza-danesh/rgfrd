//
//  space_r_bimol1_rnd.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/29/13.
//
//
 
#include "space_r_bimol1_rnd.h"

double space_r_bimol1_rnd(double r, void* params){

struct space_r_bimol1_parameters *p  = (struct space_r_bimol1_parameters *) params; // this line transform the void pointer to a type of struct displacement_irr_params
	
	double result = 0.0;
	double ka =p->ka;
	double D = p->D;
	double sigma = p->sigma;
	double r0 = p->r0;
	double rnd = p->rnd;	
	double t = p->sampled_time;

	result = space_r_bimol1(r, t, ka, D, sigma, r0) - rnd;	

	return result;
}