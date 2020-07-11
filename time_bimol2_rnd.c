//
//  time_bimol2_rnd.c
//  
//
//  Created by Zahedeh Bashardanesh on 05/26/14.
//
//
 
#include "time_bimol2_rnd.h"





double time_bimol2_rnd(double t , void* params){

	struct time_bimol2_parameters *p  = (struct time_bimol2_parameters *) params;
	
	double result;
	double rnd = p->rnd;	
	It_coeff_struct It_coeff = p->It_coeff;
	Qt_coeff_struct Qt_coeff = p->Qt_coeff;


   	result = time_bimol2(t, It_coeff, Qt_coeff) - rnd;

	
	return result;
}