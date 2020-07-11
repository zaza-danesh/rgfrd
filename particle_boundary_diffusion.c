//
//  particle_boundary_diffusion.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#include "particle_boundary_diffusion.h"

bool myfunction (particle_boundary_diffusion_time_struct struct1, particle_boundary_diffusion_time_struct struct2) { return (struct1.t_w < struct2.t_w); }

std::vector<particle_boundary_diffusion_time_struct> particle_boundary_diffusion(std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<std::vector<double> > &dij_matrix, double H, double sphere_radius, double R_curv, int k ){
	int np = particle_list.size();
	int i = -1;
	int j = -1;
	double Di = 0.0;
	double H2 = H*H;
	double shell_dist = 0.0;

	particle_boundary_diffusion_time_struct new_element;
	std::vector<particle_boundary_diffusion_time_struct> particle_boundary_diffusion_time_array;


	for (i = 0 ; i < np ; i++){
        Di = species_list[particle_list[i].index].DiffCoeff;
        shell_dist = particle_list[i].r + species_list[particle_list[i].index].radius;
        shell_dist = sphere_radius - shell_dist + sphere_radius/R_curv;
        new_element.pi = i;
        new_element.si = particle_list[i].index;
        new_element.t_w = shell_dist*shell_dist/(6*H2*Di);
        particle_list[i].t_w = new_element.t_w;
        new_element.dist = shell_dist;           
        particle_boundary_diffusion_time_array.push_back(new_element);

    }
    std::sort(particle_boundary_diffusion_time_array.begin(), particle_boundary_diffusion_time_array.end(), myfunction);
    return particle_boundary_diffusion_time_array;

}
