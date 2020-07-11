
//  check_max_time_step.h
//  
//
//  Created by Zahedeh Bashardanesh on 05/11/15.
//

//

#ifndef CHECK_MAX_TIME_STEP_H
#define CHECK_MAX_TIME_STEP_H



#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif

#include "mystructure.h"
int check_max_time_step(double& max_time_step,  int& short_time_counter);

#endif






