//
//  adjust_particles.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#ifndef ADJUST_PARTICLES_H
#define ADJUST_PARTICLES_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif

#include "mystructure.h"

void adjust_particles(std::vector<Species>& species_list, std::vector<Particle>& particle_list, double sphere_radius, int k);

#endif
