//
//  save_data.h
//  
//
//  Created by Zahedeh Bashardanesh on 04/14/15.
//
//

#ifndef SAVE_DATA_H
#define SAVE_DATA_H


#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"



void save_data( std::vector<Particle>& particle_list, std::string event, double t_previous, double t, int k);


#endif
