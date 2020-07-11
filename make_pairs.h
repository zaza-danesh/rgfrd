//
//  make_pairs.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#ifndef MAKE_PAIRS_H
#define MAKE_PAIRS_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"

std::vector<pair_struct> make_pairs(char& three_body_problem, int& number_of_pairs, double& max_diffusion_time_step, double& b_value, std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<particle_particle_diffusion_time_struct>& particle_particle_diffusion_time_array, double max_time_step_temp, double sphere_radius, double H, int k );

#endif
