//
//  get_system.c
//  
//
//  Created by Zahedeh Bashardanesh on 01/17/14.
//
//
 

#include "get_system.h"



void get_system(double& sphere_radius, int& np){

	std::string rdnt;
	std::ifstream infile;
	infile.open(DATA_PATH SYSTEM);
	infile >> rdnt;
	infile >> rdnt;
	infile >> np;
	infile >> rdnt;
	infile >> sphere_radius;
	infile.close();

}