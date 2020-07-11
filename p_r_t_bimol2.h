//
//  p_r_t_bimol2.h
//  
//
//  Created by Zahedeh Bashardanesh on 05/27/14.
//
//

#ifndef P_R_T_BIMOL2_H
#define P_R_T_BIMOL2_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"





double p_r_t_bimol2(double r, double t, double ka, double D, double sigma, double b, double kd, double gamma_b, double gamma_u, double real_root, double* char_roots );


#endif