//
//  particle_particle_diffusion.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#include "particle_particle_diffusion.h"

bool myfunction (particle_particle_diffusion_time_struct struct1, particle_particle_diffusion_time_struct struct2) { return (struct1.t_m < struct2.t_m); }


std::vector<particle_particle_diffusion_time_struct> particle_particle_diffusion(std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<std::vector<double> > &dij_matrix, double H, int k ){

	int np = particle_list.size();
	int i = -1;
	int j = -1;
	double Di = 0.0;
	double Dj = 0.0;
	double Di_p_Dj = 0.0;
	double H2 = H*H;

particle_particle_diffusion_time_struct new_element;
std::vector<particle_particle_diffusion_time_struct> particle_particle_diffusion_time_array;

	for (i = 0 ; i < np ; i++){
        Di = species_list[particle_list[i].index].DiffCoeff;
        for (j = i+1 ; j < np ; j++){
            Di_p_Dj = Di + species_list[particle_list[j].index].DiffCoeff; 
            new_element.pi_1 = i;
            new_element.pi_2 = j;
            new_element.t_m = dij_matrix[i][j]*dij_matrix[i][j]/(6*H2*Di_p_Dj);
            new_element.dist = dij_matrix[i][j];           
            particle_particle_diffusion_time_array.push_back(new_element);
        }
    }

    std::sort(particle_particle_diffusion_time_array.begin(), particle_particle_diffusion_time_array.end(), myfunction);
    return particle_particle_diffusion_time_array;
}