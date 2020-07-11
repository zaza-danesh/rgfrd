//
//  find_state.h
//  
//
//  Created by Zahedeh Bashardanesh on 11/20/13.
//
//

#ifndef FIND_STATE_H
#define FIND_STATE_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"

// #include "find_roots.h"
#include "Qt_bimol2.h"
#include "It_bimol2.h"
#include "sttoint.h"
#include "Jt_bimol2.h"


std::string find_state_intermediate(double rnd, double sampled_time, phys_parameters_bimol2 physparams, double p1, double p2);

std::string find_state_exit(double rnd, double sampled_time, phys_parameters_bimol2 physparams, double It, double Qt, double Jt);
#endif