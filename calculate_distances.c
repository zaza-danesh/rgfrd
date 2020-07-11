//
//  calculate_distances.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#include "calculate_distances.h"

char calculate_distances(std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<std::vector<double> > &dij_matrix, int np, int k){

	char negative_distances = 'F';
	// int np = particle_list.size();
	int i = -1;
	int j = -1;

	for ( i = 0 ; i < np ; i++ ){
        particle_list[i].r = sqrt(particle_list[i].x*particle_list[i].x + particle_list[i].y*particle_list[i].y + particle_list[i].z*particle_list[i].z);
        if (particle_list[i].r + species_list[particle_list[i].index].radius > 1.0 ){
        }
        particle_list[i].phi = atan(particle_list[i].y/particle_list[i].x);
        particle_list[i].theta = acos(particle_list[i].z/particle_list[i].r);
        for ( j = i+1 ; j < np ; j++ ){
            dij_matrix[i][j] = sqrt(((particle_list[i].x)- (particle_list[j].x))*((particle_list[i].x)- (particle_list[j].x)) + ((particle_list[i].y)- (particle_list[j].y))* ((particle_list[i].y)- (particle_list[j].y)) + ((particle_list[i].z)- (particle_list[j].z))*((particle_list[i].z)- (particle_list[j].z))) - (species_list[particle_list[i].index].radius + species_list[particle_list[j].index].radius);            
            if (dij_matrix[i][j] < 0 ){
                negative_distances = 'T';                
            }
        }
    }
	return negative_distances;
	
}