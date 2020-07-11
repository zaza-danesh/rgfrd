//
//  Find_Roots.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/03/11.
//
//

#ifndef FINDROOTS_H
#define FINDROOTS_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif

#include "mystructure.h"
#include "make_parameters.h"
#include "sttoint.h"
#include "find_roots.h"

std::vector<Roots> Find_Roots(std::vector<Species>& species_list, double b_value);

#endif
