//
//  Jt_coeff.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/05/06.
//
//

#ifndef JT_COEFF_H
#define JT_COEFF_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"
#include "sttoint.h"



std::vector<Jt_coeff_struct> Jt_coeff(std::vector<Species>& species_list, std::vector<Roots>& species_roots, double b_value );


#endif