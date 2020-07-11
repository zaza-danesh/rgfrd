//
//  draw_space.h
//  
//
//  Created by Zahedeh Bashardanesh on 10/30/13.
//
//

#ifndef DRAW_SPACE_H
#define DRAW_SPACE_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"

#include "survival_bimol1.h"
#include "space_r_bimol1.h"
#include "space_r_bimol1_rnd.h"
#include "space_theta_bimol1.h"
#include "space_theta_bimol1_rnd.h"
#include "f_free.h"
#include "f_corr.h"
#include "pp_free.h"
#include "pp_corr.h"
#include "survival_bimol2.h"
#include "space_r_bimol2.h"
#include "space_r_bimol2_rnd.h"
#include "find_roots.h"


double draw_space_r_bimol1(double rnd, double t, phys_parameters_bimol1 physparams);

double draw_space_theta_bimol1(double rnd, double r,  double t, phys_parameters_bimol1 physparams);
double draw_space_theta_bimol1_test(int k, double rnd, double r,  double t, phys_parameters_bimol1 physparams);

double draw_space_r_bimol2(double rnd, double t,  phys_parameters_bimol2 physparams, space_r_coeff_struct space_r_coeff,  double psurv);




#endif