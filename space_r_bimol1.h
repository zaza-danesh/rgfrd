//
//  space_r_bimol1.h
//  
//
//  Created by Zahedeh Bashardanesh on 10/29/13.
//
//

#ifndef SPACE_R_BIMOL1_H
#define SPACE_R_BIMOL1_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"



double space_r_bimol1(double r, double t , double ka, double D, double sigma, double r0);

struct space_r_bimol1_parameters {
	double ka;
	double D;
	double sigma;
	double r0;
	double rnd; 
	double sampled_time;

};

#endif 
