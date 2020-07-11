//
//  gaussian_distributing.c
//  
//
//  Created by Zahedeh Bashardanesh on 09/23/14.
//
//

 
#include "gaussian_distributing.h"

void gaussian_distributing(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi, double sampled_time, gsl_rng * gen){


	double D_i = 0.0;
	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;


// sampling the new position for the chosen particle
	D_i = species_list[particle_list[pi].index].DiffCoeff;
	dx = gsl_ran_gaussian(gen, sqrt(2*D_i*sampled_time));                
	dy = gsl_ran_gaussian(gen, sqrt(2*D_i*sampled_time));
	dz = gsl_ran_gaussian(gen, sqrt(2*D_i*sampled_time));
	particle_list[pi].x = particle_list[pi].x  + dx;
	particle_list[pi].y = particle_list[pi].y  + dy;
	particle_list[pi].z = particle_list[pi].z  + dz;

	particle_list[pi].r = sqrt(particle_list[pi].x*particle_list[pi].x + particle_list[pi].y*particle_list[pi].y + particle_list[pi].z*particle_list[pi].z);
	if(particle_list[pi].x < 0.0){
		particle_list[pi].phi = atan(particle_list[pi].y/particle_list[pi].x) + Pi;	
	}
	else{
		particle_list[pi].phi = atan(particle_list[pi].y/particle_list[pi].x);
	}
	particle_list[pi].theta = acos(particle_list[pi].z/particle_list[pi].r);
	particle_list[pi].update_flag ++;

	
}