
//  p_corr.h
//  
//
//  Created by Zahedeh Bashardanesh on 11/06/13.
//
//

#ifndef P_CORR_H
#define P_CORR_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"
#include "make_tables.h"

double p_corr(double theta, double r, double t, phys_parameters_bimol1 );





#endif