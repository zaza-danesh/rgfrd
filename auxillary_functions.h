//
//  auxillary_functions.h
//  
//
//  Created by Zahedeh Bashardanesh on 10/28/13.
//
//

#ifndef auxillary_functions_H
#define auxillary_functions_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"




double expsq_erfc(double );
double W(double, double);

double spherical_j_bessel(unsigned int, double);
double spherical_y_bessel(unsigned int, double);
double j_small_n(unsigned int, double);
double y_small_n(unsigned int, double);

#endif


