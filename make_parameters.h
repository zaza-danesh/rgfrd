//
//  make_parameters.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/25.
//
//

#ifndef MAKE_PARAMETERS_H
#define MAKE_PARAMETERS_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif

#include "mystructure.h"
phys_parameters_bimol1 make_parameters_bimol1( double ka, double D, double sigma,  double r0 );

phys_parameters_bimol2 make_parameters_bimol2(double ka, double D, double sigma, double b, double kd, double gamma_b, double gamma_u);


#endif
