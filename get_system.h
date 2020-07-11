//
//  get_system.h
//  
//
//  Created by Zahedeh Bashardanesh on 01/17/14.
//
//



#ifndef GET_SYSTEM_H
#define GET_SYSTEM_H



#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"


void get_system(double& sphere_radius, int& initial_number_of_particles);


#endif