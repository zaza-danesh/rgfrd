//  save_data.c
//  
//
//  Created by Zahedeh Bashardanesh on 04/08/15.


//  Zahedeh: This function save the particles index, the sampled time at time step k, the position of all particles available in the system and also the radius and radius^2! Further the radius multiplied by the time step!
//
 
#include "save_data.h"


void save_data(std::vector<Particle>& particle_list, std::string event, double t_previous, double t, int k){


    int np = particle_list.size();
    std::ofstream  number_of_particles_file, time_step_file;
    number_of_particles_file.open(DATA_OUTPUT_PATH NUMBER_OF_PARTICLES_FILE_NAME, std::ios_base::app);
    time_step_file.open(DATA_OUTPUT_PATH TIME_STEP_FILE_NAME, std::ios_base::app);

    number_of_particles_file << np << std::endl;
    time_step_file << t - t_previous << std::endl;
    
    time_step_file.close();
    number_of_particles_file.close();

    std::ofstream particles_n_file, particles_n2_file;
    std::ofstream particles_r_file, particles_r2_file;
    std::ofstream particles_x_file, particles_y_file, particles_z_file;

    
    particles_n_file.open(DATA_OUTPUT_PATH PARTICLES_N_FILE_NAME, std::ios_base::app);
    particles_n2_file.open(DATA_OUTPUT_PATH PARTICLES_N2_FILE_NAME, std::ios_base::app);
    particles_r_file.open(DATA_OUTPUT_PATH PARTICLES_R_FILE_NAME, std::ios_base::app);
    particles_r2_file.open(DATA_OUTPUT_PATH PARTICLES_R2_FILE_NAME, std::ios_base::app);
    particles_x_file.open(DATA_OUTPUT_PATH PARTICLES_X_FILE_NAME, std::ios_base::app);
    particles_y_file.open(DATA_OUTPUT_PATH PARTICLES_Y_FILE_NAME, std::ios_base::app);
    particles_z_file.open(DATA_OUTPUT_PATH PARTICLES_Z_FILE_NAME, std::ios_base::app);

    
    int i = 0;
    int p1_n = 0;
    int p2_n = 0;
    int p3_n = 0;
    int p4_n = 0;
    int p5_n = 0;
    int p6_n = 0;
    int flag_p1 = 0;
    int flag_p2 = 0;
    int flag_p3 = 0;
    int flag_p4 = 0;
    int flag_p5 = 0;
    int flag_p6 = 0;

    double p1_r = 0;
    double p2_r = 0;
    double p3_r = 0;
    double p4_r = 0;
    double p5_r = 0;
    double p6_r = 0;
    double delta_t = t - t_previous;


    double p1_x = 0.0;
    double p1_y = 0.0;
    double p1_z = 0.0;
    double p2_x = 0.0;
    double p2_y = 0.0;
    double p2_z = 0.0;
    double p3_x = 0.0;
    double p3_y = 0.0;
    double p3_z = 0.0;
    double p4_x = 0.0;
    double p4_y = 0.0;
    double p4_z = 0.0;
    double p5_x = 0.0;
    double p5_y = 0.0;
    double p5_z = 0.0;
    double p6_x = 0.0;
    double p6_y = 0.0;
    double p6_z = 0.0;
 

    for (i = 0 ; i < np ; i++){
        if (particle_list[i].index == 1){
            if(flag_p1 != 0){
                p1_n++;    
            }
            else{
                p1_r = particle_list[i].r;
                p1_x = particle_list[i].x;
                p1_y = particle_list[i].y;
                p1_z = particle_list[i].z;
                p1_n++;
                flag_p1 = 1;
            }     
        }
        else if (particle_list[i].index == 2){
            if (flag_p2 != 0)
            {
                p2_n++;
            }
            else{

                p2_r = particle_list[i].r;
                p2_x = particle_list[i].x;
                p2_y = particle_list[i].y;
                p2_z = particle_list[i].z;
                p2_n++;
                flag_p2 = 1;
            }
            
        }
        else if (particle_list[i].index == 3){
            if(flag_p3 != 0)
            {
                p3_n++;
            }
            else{
                p3_r = particle_list[i].r;
                p3_x = particle_list[i].x;
                p3_y = particle_list[i].y;
                p3_z = particle_list[i].z;
                p3_n++;
                flag_p3 = 1;
            }
        }
        else if (particle_list[i].index == 4){
            if (flag_p4 != 0){
                p4_n++;
            }
            else{
                p4_r = particle_list[i].r;
                p4_x = particle_list[i].x;
                p4_y = particle_list[i].y;
                p4_z = particle_list[i].z;
                p4_n++;
                flag_p4 = 1;
            }

        }
        else if (particle_list[i].index == 5){
            if (flag_p5 != 0){
                p5_n++;
            }
            else{
                p5_r = particle_list[i].r;
                p5_x = particle_list[i].x;
                p5_y = particle_list[i].y;
                p5_z = particle_list[i].z;
                p5_n++;
                flag_p5 = 1;
            }
            
        }
        else if (particle_list[i].index == 6){
            if (flag_p6 != 0){
                p6_n++;
            }
            else{
                p6_r = particle_list[i].r;
                p6_x = particle_list[i].x;
                p6_y = particle_list[i].y;
                p6_z = particle_list[i].z;
                p6_n++;
                flag_p6 = 1;
            }
            
        }
        else{
            // std::cout << __FILE__ << SP << __LINE__ << SP << "error: nonspecific species" << std::endl;
        }

    }
 


    particles_n_file << p1_n << SP << p2_n << SP << p3_n << SP << p4_n << SP << p5_n << SP << p6_n << std::endl;

    particles_x_file << p1_x << SP << p2_x << SP << p3_x << SP << p4_x << SP << p5_x <<  SP << p6_x <<  std::endl;

    particles_y_file << p1_y << SP << p2_y << SP << p3_y << SP << p4_y <<SP << p5_y <<SP << p6_y << std::endl;

    particles_z_file << p1_z << SP << p2_z << SP << p3_z << SP << p4_z << SP << p5_z <<SP << p6_z << std::endl;

    particles_n2_file << p1_n * p1_n << SP << p2_n * p2_n << SP << p3_n * p3_n << SP << p4_n * p4_n <<  SP << p5_n * p5_n <<  SP << p6_n * p6_n <<  std::endl;

    particles_r_file << p1_r << SP << p2_r << SP << p3_r << SP << p4_r << SP << p5_r <<  SP << p6_r << std::endl;

    particles_r2_file << p1_r * p1_r << SP << p2_r * p2_r << SP << p3_r * p3_r << SP << p4_r * p4_r << SP << p5_r * p5_r <<  SP << p6_r * p6_r << std::endl;


    particles_n_file.close();
    particles_n2_file.close();
    particles_r_file.close();
    particles_r2_file.close();
    particles_x_file.close();
    particles_y_file.close();
    particles_z_file.close();



}    
