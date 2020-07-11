//
//  space_r_bimol1_rnd.h
//  
//
//  Created by Zahedeh Bashardanesh on 10/29/13.
//
//


#ifndef SPACE_R_BIMOL1_RND_H
#define SPACE_R_BIMOL1_RND_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"
#include  "space_r_bimol1.h"

double space_r_bimol1_rnd(double r, void* params);

#endif
