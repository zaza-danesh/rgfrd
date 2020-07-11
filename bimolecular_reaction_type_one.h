//  bimolecular_reaction_type_one.h
//  
//
//  Created by Zahedeh Bashardanesh on 11/07/14.
//
//

#ifndef BIMOLECULAR_REACTION_TYPE_ONE_H
#define BIMOLECULAR_REACTION_TYPE_ONE_H



#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "sttoint.h"

void bimol1_bound(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int p1, int p2, double sampled_time, gsl_rng * gen);

void bimol1_unbound(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int p1, int p2, double sampled_time, double sampled_r, double sampled_theta, gsl_rng * gen);


#endif






