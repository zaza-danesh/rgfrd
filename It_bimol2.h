//
//  It_bimol2.h
//  
//
//  Created by Zahedeh Bashardanesh on 05/25/14.
//
//

#ifndef IT_BIMOL2_H
#define IT_BIMOL2_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"


double It_bimol2(double t, It_coeff_struct It_coeff );


#endif