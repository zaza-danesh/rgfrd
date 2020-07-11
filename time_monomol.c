
//
//  time_monomol_distribution.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/05/13.
//
//
 
#include "time_monomol.h"

double time_monomol_distribution(double t, double kd){

	double result = 0.0;

	result = 1-exp(-kd*t);
	
	return result;

}

