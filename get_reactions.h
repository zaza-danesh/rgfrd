//
//  get_reactions.h
//  
//
//  Created by Zahedeh Bashardanesh on 11/19/15.
//
//

#ifndef GET_REACTIONS_H
#define GET_REACTIONS_H



#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif
#include "mystructure.h"

void get_reactions(int ns, std::vector<Species>& species_list, double ka1, double kd1, double k1, double ka2, double kd2, double k2);


#endif