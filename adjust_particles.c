//
//  adjust_particles.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#include "adjust_particles.h"


void adjust_particles(std::vector<Species>& species_list, std::vector<Particle>& particle_list, double sphere_radius, int k){

    int np = particle_list.size();
    int i = -1;
    int j = -1;
    int kk = -1;
    double d_ij = 0.0;
    double d_ij2 = 0.0;
    double sigma = 0.0;
    double sigma2 = 0.0;
    double den = 0.0;
    double deltaX = 0.0;
    double deltaY = 0.0;
    double deltaZ = 0.0;
    double epsilon_m = 2e-4;
    double epsilon_w = 1e-4;
    double force = 0.0;

    for (kk = 0 ; kk < ADJUSTMENT_STEPS_NUMBER ; kk++){
        for(i = 0 ; i < np ; i++){
            particle_list[i].deltaX = 0.0;
            particle_list[i].deltaY = 0.0;              
            particle_list[i].deltaZ = 0.0;
        }
        for (i = 0 ; i < np ; i++){
            for (j = i+1 ; j < np ; j++){
                d_ij = sqrt((particle_list[i].x - particle_list[j].x)*(particle_list[i].x - particle_list[j].x) + (particle_list[i].y - particle_list[j].y)*(particle_list[i].y - particle_list[j].y) + (particle_list[i].z - particle_list[j].z)*(particle_list[i].z - particle_list[j].z))- (species_list[particle_list[i].index].radius + species_list[particle_list[j].index].radius);
                d_ij2 = d_ij*d_ij;
                sigma = species_list[particle_list[i].index].radius + species_list[particle_list[j].index].radius;
                sigma2 = sigma*sigma;
                //x-direction
                deltaX = sigma2*(particle_list[i].x - particle_list[j].x)/(sigma2 + d_ij2);
                particle_list[i].deltaX = particle_list[i].deltaX + deltaX;
                particle_list[j].deltaX = particle_list[j].deltaX - deltaX;
                //y-direction
                deltaY = sigma2*(particle_list[i].y - particle_list[j].y)/(sigma2 + d_ij2);
                particle_list[i].deltaY = particle_list[i].deltaY + deltaY;
                particle_list[j].deltaY = particle_list[j].deltaY - deltaY;
                //z-direction
                deltaZ = sigma2*(particle_list[i].z - particle_list[j].z)/(sigma2 + d_ij2);
                particle_list[i].deltaZ = particle_list[i].deltaZ + deltaZ;
                particle_list[j].deltaZ = particle_list[j].deltaZ - deltaZ;
            }
        }
        for (i = 0 ; i < np ; i++){
            epsilon_w = 1e-4;
            particle_list[i].r = sqrt(particle_list[i].x * particle_list[i].x + particle_list[i].y*particle_list[i].y +  particle_list[i].z*particle_list[i].z );
            force = species_list[particle_list[i].index].radius *  species_list[particle_list[i].index].radius/((sphere_radius - species_list[particle_list[i].index].radius) - particle_list[i].r );
            if(fabs(force)*epsilon_w > species_list[particle_list[i].index].radius ){
                epsilon_w = epsilon_w/ADJUSTMENT_STEPS_NUMBER;
                }
            particle_list[i].x = particle_list[i].x + epsilon_m*particle_list[i].deltaX - epsilon_w*force*particle_list[i].x/particle_list[i].r;
            particle_list[i].y = particle_list[i].y + epsilon_m*particle_list[i].deltaY - epsilon_w*force*particle_list[i].y/particle_list[i].r;
            particle_list[i].z = particle_list[i].z + epsilon_m*particle_list[i].deltaZ - epsilon_w*force*particle_list[i].z/particle_list[i].r;
        }
    }

    
}

