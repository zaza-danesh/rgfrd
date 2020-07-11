//
//  draw_event_time.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#ifndef DRAW_EVENT_TIME_H
#define DRAW_EVENT_TIME_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif


#include "mystructure.h"
#include "get_parameters.h"
#include "make_parameters.h"
#include "draw_time.h"

event_struct draw_event_time( std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<pair_struct>& pair_list, int number_of_pairs, std::vector<It_coeff_struct> It_coeff_arr, std::vector<Qt_coeff_struct> Qt_coeff_arr, double b_value, gsl_rng * gen, int k);

#endif
