//
//  f_free.h
//  
//
//  Created by Zahedeh Bashardanesh on 11/06/13.
//
//

#ifndef F_FREE_H
#define F_FREE_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"

// double integrand_J0(double , double , double , phys_parameters_association );

// double integrate_J0(double , void* params);

// double J0( double , double , gsl_integration_workspace* , double, phys_parameters_association );

double f_free(double r, double t, phys_parameters_bimol1 physparams);


// struct integrand_J0_parameters { 
// 	double r;
// 	double t;
// 	double ka;
// 	double D;
// 	double sigma; 
// 	double r0;
// };


#endif	