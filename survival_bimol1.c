//
//  survival_bimol1.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/29/13.
//
//
 
#include "survival_bimol1.h"

double survival_bimol1(double t, double ka, double D, double sigma, double r0){

	
	return (1 - time_bimol1(t,  ka, D, sigma, r0));
}

