//
//  time_monomol_rnd.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/05/13.
//
//
 
#include "time_monomol_rnd.h"



double time_monomol_rnd(double t, void* params){

struct time_monomol_params *p  = (struct time_monomol_params *) params; // this line transform the void pointer to a type of struct time_irr_params

	double result = 0.0;

	double kd = p->kd;	
	double rnd = p->rnd;

	
	result = time_monomol_distribution( t,  kd ) - rnd;


	return result;
}

