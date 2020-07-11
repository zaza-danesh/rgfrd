
//  reflection.c
//  
//
//  Created by Zahedeh Bashardanesh on 04/08/15.


//  Zahedeh: This function checks if the particles are all inside the container, if not it reflects it back inside perpendicularly! 
//
 
#include "reflection.h"


void reflection(std::vector<Species>& species_list, std::vector<Particle>& particle_list, double sphere_radius, double& ref_flag, double final_time, int k){


    int i = 0;
    int reflection_flag = 0;
    double r_temp = 0.0;
    double shell_r = 0.0;
    std::string event = "--";

    for (i = 0 ; i < particle_list.size() ; i++){
        r_temp = particle_list[i].r;
        shell_r = particle_list[i].r +  species_list[particle_list[i].index].radius;
        if (shell_r>sphere_radius ){
            if (shell_r > 2*sphere_radius){
                ref_flag = -1;
            }
            else{//reflecting the particle back
                reflection_flag ++;
                particle_list[i].reflect_flag = 1;
                shell_r = 2*sphere_radius - shell_r;
                particle_list[i].r = shell_r - species_list[particle_list[i].index].radius;
                particle_list[i].x = particle_list[i].x*particle_list[i].r/r_temp;
                particle_list[i].y = particle_list[i].y*particle_list[i].r/r_temp;
                particle_list[i].z = particle_list[i].z*particle_list[i].r/r_temp;
                particle_list[i].phi = atan(particle_list[i].y/particle_list[i].x);
                particle_list[i].theta = acos(particle_list[i].z/particle_list[i].r);
            }
        }
    }
}


