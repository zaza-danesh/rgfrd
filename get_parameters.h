//
//  get_parameters.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/25.
//
//

#ifndef GET_PARAMETERS_H
#define GET_PARAMETERS_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif

#include "mystructure.h"
#include "sttoint.h"

void get_parameters_bimol1(std::vector<Species>& species_list, int I, int J, double* ka, double* D, double* sigma);

void get_parameters_bimol2(std::vector<Species>& species_list, int I, int J, double* ka, double* D, double* sigma, double* kd, double* gamma_b, double* gamma_u);


#endif
