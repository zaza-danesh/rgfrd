//
//  get_species.h
//  
//
//  Created by Zahedeh Bashardanesh on 01/17/14.
//
//



#ifndef GET_SPECIES_H
#define GET_SPECIES_H



#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"


void get_species(int number_of_species, Species* species_list);
int get_number_of_species();



#endif