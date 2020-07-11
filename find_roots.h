//
//  find_roots.h
//  
//
//  Created by Zahedeh Bashardanesh on 12/13/13
//  Modified by Zahedeh Bashardanesh on 11/04/14
//
//

#ifndef FIND_ROOTS_H
#define FIND_ROOTS_H


#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"




void find_roots( double& real_root, double* char_roots,  phys_parameters_bimol2 physparams);

double find_roots_hx(phys_parameters_bimol2 physparams, double x_low, double x_high);

double hx(double l, void *params);

struct hx_params {

	double kappa;
	double kappa1_r;
	double kd;
	double N_disso;
	double gamma_b;
	double gamma_u;
};

double find_roots_hy(phys_parameters_bimol2 physparams, double y_low, double y_high);

double find_roots_htan(phys_parameters_bimol2 physparams, double y_low, double y_high);

double hy(double q, void *params);

double htan(double q, void *params);


struct hy_params {

	double kappa;
	double kappa1_r;
	double kd;
	double N_disso;
	double gamma_b;
	double gamma_u;
};

double sketch_hx(double l, double kappa, double kappa1_r, double kd, double N_disso, double gamma_b, double gamma_u);

double sketch_hy(double y, double kappa, double kappa1_r, double kd, double N_disso, double gamma_b, double gamma_u);


#endif