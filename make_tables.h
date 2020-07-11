//
//  make_tables.h
//  
//
//  Created by Zahedeh Bashardanesh on 10/31/13.
//
//

#ifndef MAKE_TABLES_H
#define MAKE_TABLES_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"


// double p_theta_free(double, double, double, phys_parameters_bimol1);


// double ip_theta_free(double, double, double, phys_parameters_bimol1); 

// double p_free_max(double , double, phys_parameters_bimol1);



void makeRnTable(doubleVector& RnTable, double r, double t, phys_parameters_bimol1 physparams);

void makeLegendreTable(doubleVector const& RnTable, double theta );

double integrand_Rn(unsigned int n, double u, double r, double t, phys_parameters_bimol1 physparams);

double integrate_Rn(double u, void* params);

double Rn(unsigned int n, double r, double t, gsl_integration_workspace* workspace, double tolerance, phys_parameters_bimol1 physparams );

double p_irr(double r, double t, phys_parameters_bimol1 physparams);

double ip_theta_free(double theta,  double r, double t, phys_parameters_bimol1 physparams);

double p_free_max(double r , double t, phys_parameters_bimol1 physparams);


struct integrand_Rn_parameters { 
	unsigned int n;
	double r;
	double t;
	double ka;
	double D;
	double sigma; 
	double r0;
};

#endif