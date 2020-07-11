//
//  initial_distribution.c
//  
//
//  Created by Zahedeh Bashardanesh on 01/13/14.
//
//
 


#include "initial_distribution.h"


std::vector<Particle> initial_distribution(std::vector<Species>& species_list,  std::vector<std::vector<double> > &dij_matrix, int np,  double sphere_radius, gsl_rng * gen, int k){


    char negative_distance = 'T';
    int fin = species_list[1].quantity;

    int init = 0;   
    int i = 0;
    int ii=0;
    int j = 0;
    int kk = 0;
    int ns = species_list.size();
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double v = 0.0;
    double r = 0.0;
    double d_ij = 0.0;
    double sphere_factor = 2.0/3.0;
    Particle new_particle;
    std::vector<Particle> particle_list(np);

    while (negative_distance == 'T'){
        init = 0;
        fin = species_list[1].quantity; 
        kk++;
        for (i = 0 ; i < np ; i++){
            r = sphere_factor*sphere_radius *pow(gsl_rng_uniform(gen), 1.0/3.0);
            x = gsl_ran_gaussian(gen, 1.0);
            y = gsl_ran_gaussian(gen, 1.0);
            z = gsl_ran_gaussian(gen, 1.0);
            v = sqrt(x*x + y*y + z*z);
            particle_list[i].x = r*x/v;
            particle_list[i].y = r*y/v;
            particle_list[i].z = r*z/v;
            particle_list[i].r = sqrt(particle_list[i].x*particle_list[i].x + particle_list[i].y*particle_list[i].y + particle_list[i].z*particle_list[i].z);
            particle_list[i].phi = atan(particle_list[i].y/particle_list[i].x);
            particle_list[i].theta = acos(particle_list[i].z/particle_list[i].r);
        }//end of distributing particles

        for (j = 1 ; j < ns ; j++){
           for (ii = init ; ii < fin ; ii++ ){
            // std::cout << "ii: " << ii << std::endl;
                particle_list[ii].index = j;

            }
            init = fin;
            fin = fin + species_list[j+1].quantity;
        }
        negative_distance = calculate_distances( species_list, particle_list, dij_matrix, np, k);
    }

    return particle_list;

}