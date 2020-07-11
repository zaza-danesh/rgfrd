//
//  Jt_bimol2.h
//  
//
//  Created by Zahedeh Bashardanesh on 02/25/15.
//
//

#ifndef JT_BIMOL2_H
#define JT_BIMOL2_H

// #include "/Users/zahedeh/MyInclude/mystlibgpp.h"
#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"


double Jt_bimol2(double t,  Jt_coeff_struct Jt_coeff);


#endif