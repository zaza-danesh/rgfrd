//
//  draw_time.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#ifndef DRAW_TIME_H
#define DRAW_TIME_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif

#include "mystructure.h"
#include "time_bimol1_rnd.h"
#include "time_bimol2_rnd.h"
#include "time_monomol_rnd.h"
#include "find_roots.h"

double draw_time_bimol1(double rnd,  phys_parameters_bimol1 physparams);

double draw_time_bimol2(int k, double rnd,  phys_parameters_bimol2 physparams, It_coeff_struct It_coeff, Qt_coeff_struct Qt_coeff);

double draw_time_monomol(double rnd, double physparams_exp);

#endif
