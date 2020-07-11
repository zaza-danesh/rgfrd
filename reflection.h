//
//  reflection.h
//  
//
//  Created by Zahedeh Bashardanesh on 04/08/15.
//
//

#ifndef REFLECTION_H
#define REFLECTION_H


#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"



void reflection(std::vector<Species>& species_list, std::vector<Particle>& particle_list, double sphere_radius, double& ref_flag, double final_time, int k);


#endif
