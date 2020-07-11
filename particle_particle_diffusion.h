//
//  particle_particle_diffusion.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#ifndef PARTICLE_PARTICLE_DIFFUSION_H
#define PARTICLE_PARTICLE_DIFFUSION_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif

#include "mystructure.h"


std::vector<particle_particle_diffusion_time_struct> particle_particle_diffusion(std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<std::vector<double> > &dij_matrix, double H, int k );

#endif
