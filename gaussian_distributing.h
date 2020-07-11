//
//  gaussian_distributing.h
//  
//
//  Created by Zahedeh Bashardanesh on 09/23/14.
//
//

#ifndef GAUSSIAN_DISTRIBUTING_H
#define GAUSSIAN_DISTRIBUTING_H



#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "sttoint.h"
#include "draw_time.h"

void gaussian_distributing( std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int I_transf, double sampled_time, gsl_rng * gen);

#endif