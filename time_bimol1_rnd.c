//
//  time_bimol1_rnd.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/28/13.
//
//
 

#include "time_bimol1_rnd.h"



double time_bimol1_rnd(double t, void* params){

struct time_bimol1_parameters *p  = (struct time_bimol1_parameters *) params; // this line transform the void pointer to a type of struct time_bimol1_parameters

	double result = 0.0;
	double ka =p->ka;
	double D = p->D;
	double sigma = p->sigma;
	double r0 = p->r0;
	double rnd = p->rnd;
	

 

	result = time_bimol1(t, ka, D, sigma, r0) - rnd;


	return result;
}

