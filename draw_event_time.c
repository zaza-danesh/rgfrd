//
//  draw_event_time.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#include "draw_event_time.h"

event_struct draw_event_time( std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<pair_struct>& pair_list, int number_of_pairs, std::vector<It_coeff_struct> It_coeff_arr, std::vector<Qt_coeff_struct> Qt_coeff_arr, double b_value, gsl_rng * gen, int k){

    int i = -1;
    int ii = -1;
    int index = -1;
    int np = particle_list.size();
    int pair_index = -1;
    int pi_bimol2 = -1;
    int pi_monomol = -1;
    double rnd_bimol1 = 0.0;
    double rnd_bimol2 = 0.0;
    double rnd_monomol = 0.0;
    double t_bimol1 = 0.0;
    double temp_t_bimol1 = INF;
    double t_bimol2 = 0.0;
    double temp_t_bimol2 = INF;
    double t_monomol = 0.0;
    double temp_t_monomol = INF;
    double ka = 0.0;
    double D = 0.0;
    double sigma = 0.0;
    double r0 = 0.0;
    double factor = 0.0;
    double b = 0.0;
    double kd = 0.0;
    double gamma_b = 0.0;
    double gamma_u = 0.0;
    double physparams_monomol;
    double temp_physparams_monomol;
    event_struct reaction_event;
    phys_parameters_bimol1 physparams_bimol1;
    phys_parameters_bimol2 physparams_bimol2;
    phys_parameters_bimol1 temp_physparams_bimol1;
    phys_parameters_bimol2 temp_physparams_bimol2;


    for (i = 0 ; i < number_of_pairs ; i++){        
        get_parameters_bimol1(species_list, particle_list[pair_list[i].p1].index, particle_list[pair_list[i].p2].index, &ka, &D, &sigma);
        r0 = pair_list[i].dist + species_list[particle_list[pair_list[i].p1].index].radius + species_list[particle_list[pair_list[i].p2].index].radius;
        factor = sigma*ka/(r0*(ka + 4*Pi*sigma*D));
        rnd_bimol1 = gsl_rng_uniform(gen);
        if (rnd_bimol1 > factor){
            t_bimol1 = INF-(double)i;
        }
        else{
            physparams_bimol1 =  make_parameters_bimol1(ka, D, sigma, r0);
            t_bimol1 = draw_time_bimol1(rnd_bimol1, physparams_bimol1); 
        }
        if (t_bimol1 < temp_t_bimol1){
            temp_t_bimol1 = t_bimol1;
            pair_index = i;
        }
    }
    // temp_t_bimol2 = temp_t_bimol1; 
    // temp_t_monomol = temp_t_bimol1;
    for ( i = 0 ; i < np ; i++){
        if (species_list[particle_list[i].index].disso_flag != 0){ // should check for the cases that disso_flag == 1
            if (particle_list[i].pair_flag != 1){
                get_parameters_bimol2(species_list, particle_list[i].index, 0, &ka, &D, &sigma, &kd, &gamma_b, &gamma_u);

                rnd_bimol2 = gsl_rng_uniform(gen);
                if (ka != 0.0 ){
                    // b = std::min( ALPHA*sigma, b_value + sigma);
                    b = ALPHA*sigma;
                    physparams_bimol2 = make_parameters_bimol2(ka, D, sigma, b, kd, gamma_b, gamma_u);
                    for ( ii = 0 ; ii < It_coeff_arr.size() ; ii++){
                        if (It_coeff_arr[ii].index == particle_list[i].index){                            
                            index = ii;
                        }
                    }
                    t_bimol2 = draw_time_bimol2(k, rnd_bimol2, physparams_bimol2, It_coeff_arr[index],  Qt_coeff_arr[index]);
                }
                else{
                    t_bimol2 = draw_time_monomol(rnd_bimol2, kd + gamma_b);
                    b = sigma;
                    physparams_bimol2 = make_parameters_bimol2(ka, D, sigma, b, kd, gamma_b, gamma_u);
                }
                if (t_bimol2 < temp_t_bimol2){
                    temp_t_bimol2 = t_bimol2;
                    pi_bimol2 = i;
                    // temp_physparams_bimol2 = physparams_bimol2;
                }               
            }
            else if((species_list[particle_list[i].index].trans_bimol_flag != 0 || species_list[particle_list[i].index].trans_monomol_flag != 0 ) && species_list[particle_list[i].index].disso_flag == 0 ){
                rnd_monomol = gsl_rng_uniform(gen);
                t_monomol = draw_time_monomol(rnd_monomol, species_list[particle_list[i].index].total_transform_rate);
                if (t_monomol < temp_t_monomol){
                    // temp_physparams_monomol = species_list[particle_list[i].index].total_transform_rate;
                    temp_t_monomol = t_monomol;                    
                    pi_monomol = i;
                }     
            }
        }
        else{ 
            // std::cout << __FILE__ << SP << __LINE__ <<  ": Error: there is a particle with specification beyond the algorithm is developed for." << std::endl;
        }
    }

    if ( pi_monomol == -1 &&  pair_index == -1 && pi_bimol2 == -1){
        return reaction_event;
    }
    else{
        if (temp_t_bimol1 < temp_t_bimol2){
            if (temp_t_bimol1 < temp_t_monomol){
                reaction_event.event = "bimol1";
                reaction_event.t_r = temp_t_bimol1;
                reaction_event.pair_index = pair_index;
                // reaction_event.physparams_bimol1 = temp_physparams_bimol1;
            }
            else{
                reaction_event.event = "monomol";
                reaction_event.t_r = temp_t_monomol;
                reaction_event.pi_monomol = pi_monomol;
                // reaction_event.physparams_monomol = temp_physparams_monomol;
            }
        }
        else{
            if (temp_t_bimol2 < temp_t_monomol){
                reaction_event.event = "bimol2";
                reaction_event.t_r = temp_t_bimol2;
                reaction_event.pi_bimol2 = pi_bimol2;
                // reaction_event.physparams_bimol2 = temp_physparams_bimol2;
            }
            else{
                reaction_event.event = "monomol";
                reaction_event.t_r = temp_t_monomol;
                reaction_event.pi_monomol = pi_monomol;
                // reaction_event.physparams_monomol = temp_physparams_monomol;
            }
        }
        return reaction_event;
    }
    

    
}