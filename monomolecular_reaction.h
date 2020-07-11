//
//  monomolecular_reaction.h
//  
//
//  Created by Zahedeh Bashardanesh on 09/23/14.
//
//

#ifndef MONOMOLECULAR_REACTION_H
#define MONOMOLECULAR_REACTION_H



#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "sttoint.h"
#include "gaussian_distributing.h"


void monomolecular_reaction( std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_monomol, double t_sampled, gsl_rng * gen);
void monomolecular_reaction_bimol_or_monomol_product(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_monomol, gsl_rng * gen);
void monomolecular_reaction_bimol_product(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_monomol, int transfer_to_index, gsl_rng * gen);
#endif