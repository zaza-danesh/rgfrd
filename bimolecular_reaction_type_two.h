//
//  bimolecular_reaction_type_two.h
//  
//
//  Created by Zahedeh Bashardanesh on 09/05/14.
//
//

#ifndef BIMOLECULAR_REACTION_TYPE_TWO_H
#define BIMOLECULAR_REACTION_TYPE_TWO_H



#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "sttoint.h"
#include "monomolecular_reaction.h"
#include "gaussian_distributing.h"


void intermediate_unbound( std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_bimol2, int J_disso, double t_sampled, double r_sampled, gsl_rng * gen);

void separation( std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_bimol2, int J_disso, double t_sampled, double sep_dist, gsl_rng * gen);

void exit_bound(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_bimol2, double sampled_time, gsl_rng * gen);

void exit_unbound(std::vector<Species>& species_list, std::vector<Particle>& particle_list, int pi_bimol2, int J_disso, double sampled_time, double sampled_r, gsl_rng * gen);


#endif