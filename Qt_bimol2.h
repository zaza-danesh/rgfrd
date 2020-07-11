//
//  Qt_bimol2.h
//  
//
//  Created by Zahedeh Bashardanesh on 05/26/14.
//
//

#ifndef QT_BIMOL2_H
#define QT_BIMOL2_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "auxillary_functions.h"

double Qt_bimol2(double t,  Qt_coeff_struct Qt_coeff);


#endif
