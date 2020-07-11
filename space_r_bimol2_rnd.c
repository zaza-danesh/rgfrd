//
//  space_r_bimol2_rnd.c
//  
//
//  Created by Zahedeh Bashardanesh on 05/28/14.
//
//
 
#include "space_r_bimol2_rnd.h"



double space_r_bimol2_rnd(double r , void* params){

	struct space_bimol2_parameters *p  = (struct space_bimol2_parameters *) params;
	
	int i = 0;
	double result;
	double sampled_time = p->sampled_time;
	double rnd = p->rnd;
	space_r_coeff_struct space_r_coeff = p->space_r_coeff;
    
	result = space_r_bimol2(r, sampled_time, space_r_coeff) - rnd;
	
	return result;
}





