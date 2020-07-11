
//
//  initial_distribution.c
//  
//
//  Created by Zahedeh Bashardanesh on 01/13/14.
//
//


#ifndef INITAL_DISTRIBUTION_H
#define INITAL_DISTRIBUTION_H


#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"
#include "calculate_distances.h"
//#include  "getSystem.h"


std::vector<Particle> initial_distribution( std::vector<Species>& species_list,  std::vector<std::vector<double> > &dij_matrix, int ns, double sphere_radius,  gsl_rng * gen, int k);



#endif