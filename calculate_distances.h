
//  calculate_distances.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#ifndef CALCULATE_DISTANCES_H
#define CALCULATE_DISTANCES_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif

#include "mystructure.h"
char calculate_distances(std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<std::vector<double> >& dij_matrix, int np, int k);

#endif
