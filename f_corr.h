
//  f_corr.h
//  
//
//  Created by Zahedeh Bashardanesh on 11/06/13.
//
//

#ifndef F_CORR_H
#define F_CORR_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"

double f_corr(double r, double t, phys_parameters_bimol1 physparams);

double integrand_R0(double u, double r, double t, phys_parameters_bimol1 );

double integrate_R0(double u, void* params);
 
double R0( double r, double t, gsl_integration_workspace* workspace, double tolerance, phys_parameters_bimol1 physparams );




struct integrand_R0_parameters { 
	double r;
	double t;
	double ka;
	double D;
	double sigma; 
	double r0;
};



#endif