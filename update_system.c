//
//  update_system.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#include "update_system.h"


void update_system_bimolecular_reaction_type_one(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, std::vector<pair_struct>& pair_list, int pair_index,  std::vector<It_coeff_struct> It_coeff_arr, std::vector<Qt_coeff_struct> Qt_coeff_arr, std::vector<space_r_coeff_struct> space_r_coeff_arr, double b_value, double sampled_time, gsl_rng * gen){

    int i = -1;
    int ii = -1;
    int index = -1;
    int np = particle_list.size();
    int number_of_pairs = pair_list.size();
    std::string bimol2_state = "--";
    double ka = 0.0;
    double D = 0.0;
    double sigma = 0.0;
    double r0 = 0.0;
    double b = 0.0;
    double kd = 0.0;
    double gamma_b = 0.0;
    double gamma_u = 0.0;
    double rnd_state  = 0.0;
    double rnd_space = 0.0;
    double sampled_r = 0.0;
    double sampled_theta = 0.0;
    double It =  0.0;
    double Qt =  0.0;

    phys_parameters_bimol1 physparams_bimol1;
    phys_parameters_bimol2 physparams_bimol2;
    // std::ofstream  b_value_file;
    // b_value_file.open(DATA_OUTPUT_PATH B_VALUE_FILE_NAME, std::ios_base::app);
    

    bimol1_bound(species_list, particle_list, pair_list[pair_index].p1, pair_list[pair_index].p2, sampled_time, gen);

    for (i = 0 ; i < np ; i++ ){
        if (species_list[particle_list[i].index].disso_flag == 1 && particle_list[i].pair_flag != 1){
            get_parameters_bimol2(species_list, particle_list[i].index, 0, &ka, &D, &sigma, &kd, &gamma_b, &gamma_u);
            if (ka != 0.0 ){
                    // b = std::min( ALPHA*sigma, b_value + sigma);
                    b = ALPHA*sigma;
                    // b_value_file << b << std::endl;
                    physparams_bimol2 = make_parameters_bimol2(ka, D, sigma, b, kd, gamma_b, gamma_u);
                    rnd_state = gsl_rng_uniform(gen);
                    for ( ii = 0 ; ii < It_coeff_arr.size() ; ii++){
                        if (It_coeff_arr[ii].index == particle_list[i].index){
                            index = ii;
                        }
                    }
                    It =  It_bimol2(sampled_time, It_coeff_arr[index]);
                    Qt =  Qt_bimol2(sampled_time, Qt_coeff_arr[index]);

                    bimol2_state = find_state_intermediate(rnd_state, sampled_time, physparams_bimol2, It, Qt);
                if (bimol2_state == "bound"){
                    gaussian_distributing( species_list,  particle_list, i, sampled_time, gen);
                }
                else{
                    rnd_space = gsl_rng_uniform(gen); 
                    // if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
                    sampled_r = draw_space_r_bimol2(rnd_space, sampled_time, physparams_bimol2, space_r_coeff_arr[index], It);
                    intermediate_unbound(species_list,  particle_list, i, 0, sampled_time, sampled_r, gen);
                }
            }
            else{ //ka == 0
                std::cout << __FILE__ << SP << __LINE__ << ": Error: there is a particle with specification beyond the algorithm is developed for." << std::endl;
                std::cout << __FILE__ << SP << __LINE__ << ": Error: this particle should have been treated as a bimolecular transformation." << std::endl;
            }
            
        }
        if (species_list[particle_list[i].index].disso_flag != 1 && particle_list[i].pair_flag != 1 && particle_list[i].update_flag != 1) {
            gaussian_distributing( species_list,  particle_list, i, sampled_time, gen);
        }
    }
    for (i = 0 ; i < pair_index ; i++){
        rnd_space = gsl_rng_uniform(gen);
        // if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        get_parameters_bimol1(species_list, particle_list[pair_list[i].p1].index, particle_list[pair_list[i].p2].index, &ka, &D, &sigma);
        r0 = pair_list[i].dist +species_list[particle_list[pair_list[i].p1].index].radius + species_list[particle_list[pair_list[i].p2].index].radius ;
        physparams_bimol1 =  make_parameters_bimol1(ka, D, sigma, r0);
        sampled_r = draw_space_r_bimol1(rnd_space, sampled_time, physparams_bimol1);
        rnd_space = gsl_rng_uniform(gen);
        // if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        sampled_theta = draw_space_theta_bimol1(rnd_space, sampled_r, sampled_time, physparams_bimol1);
        bimol1_unbound(species_list, particle_list, pair_list[i].p1, pair_list[i].p2, sampled_time, sampled_r, sampled_theta, gen);

    }
    for (i = pair_index+1 ; i < number_of_pairs ; i++){
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        get_parameters_bimol1(species_list, particle_list[pair_list[i].p1].index, particle_list[pair_list[i].p2].index, &ka, &D, &sigma);
        r0 = pair_list[i].dist + species_list[particle_list[pair_list[i].p1].index].radius + species_list[particle_list[pair_list[i].p2].index].radius ;
        physparams_bimol1 =  make_parameters_bimol1(ka, D, sigma, r0);
        sampled_r = draw_space_r_bimol1(rnd_space, sampled_time, physparams_bimol1);
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        sampled_theta = draw_space_theta_bimol1(rnd_space, sampled_r, sampled_time, physparams_bimol1);
        bimol1_unbound(species_list, particle_list, pair_list[i].p1, pair_list[i].p2, sampled_time, sampled_r, sampled_theta, gen);
    }

    for (i = 0 ; i < np ; i++){
        if (particle_list[i].exist_flag == 0){
            particle_list.erase(particle_list.begin()+i);
            np--;
        }
    }
// b_value_file.close();

}

void update_system_bimolecular_reaction_type_two(int k, std::vector<Species>& species_list,  std::vector<Particle>& particle_list, std::vector<pair_struct>& pair_list, int pi_bimol2, std::vector<It_coeff_struct> It_coeff_arr, std::vector<Qt_coeff_struct> Qt_coeff_arr,std::vector<Jt_coeff_struct> Jt_coeff_arr, std::vector<space_r_coeff_struct> space_r_coeff_arr, double b_value, double sampled_time, gsl_rng * gen){

    int i = -1;
    int ii = -1;
    int index = -1;
    int np = particle_list.size();
    int number_of_pairs = pair_list.size();
    std::string bimol2_state_intermediate = "--";
    std::string bimol2_state_exit = "--";
    double ka = 0.0;
    double D = 0.0;
    double sigma = 0.0;
    double r0 = 0.0;
    double b = 0.0;
    double kd = 0.0;
    double gamma_b = 0.0;
    double gamma_u = 0.0;
    double rnd_state  = 0.0;
    double rnd_space = 0.0;
    double sampled_r = 0.0;
    double sampled_theta = 0.0;
    double It =  0.0;
    double Qt =  0.0;
    double Jt =  0.0;

    // std::ofstream  b_value_file;
    // b_value_file.open(DATA_OUTPUT_PATH B_VALUE_FILE_NAME, std::ios_base::app);

    phys_parameters_bimol1 physparams_bimol1;
    phys_parameters_bimol2 physparams_bimol2;
    get_parameters_bimol2(species_list, particle_list[pi_bimol2].index, 0, &ka, &D, &sigma, &kd, &gamma_b, &gamma_u);    
    // b =  std::min( ALPHA*sigma, b_value + sigma);
    b = ALPHA*sigma;
    for ( ii = 0 ; ii < It_coeff_arr.size() ; ii++){
        if (It_coeff_arr[ii].index == particle_list[pi_bimol2].index){
            index = ii;
        }
    }    
    rnd_state = gsl_rng_uniform(gen);    
    physparams_bimol2 = make_parameters_bimol2(ka, D, sigma, b, kd, gamma_b, gamma_u);
    // std::cout << "Index: " << index << std::endl;
    It =  It_bimol2(sampled_time, It_coeff_arr[index]);
    Qt =  Qt_bimol2(sampled_time, Qt_coeff_arr[index]);
    Jt =  Jt_bimol2(sampled_time, Jt_coeff_arr[index]);
    // std::cout << ka << SP << D << SP << sigma << SP << b << SP << kd<< SP << gamma_b <<SP << gamma_u << std::endl;
    // std::cout << "sampled_time: " << sampled_time << SP <<  "Jt: " << Jt << std::endl;    
    bimol2_state_exit = find_state_exit(rnd_state, sampled_time, physparams_bimol2, It, Qt, Jt);    
    if ( bimol2_state_exit == "separation" ){
        separation(species_list,  particle_list, pi_bimol2, 0, sampled_time, b, gen);        
    }
    else if (bimol2_state_exit == "bound" ){        
        exit_bound(species_list,  particle_list, pi_bimol2, sampled_time, gen);         
    }
    else if (bimol2_state_exit == "unbound" ){        
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        sampled_r = draw_space_r_bimol2(rnd_space, sampled_time, physparams_bimol2, space_r_coeff_arr[index], It);        
        exit_unbound(species_list,  particle_list, pi_bimol2, 0, sampled_time, sampled_r, gen); //to which particle it transforms/decays        
    }
    else {        
        std::cout << __FILE__ << SP << __LINE__ << ": Error: no exit event has been chosen for the particle." << std::endl;
    }

    for (i = 0 ; i < np ; i++ ){        
        if (species_list[particle_list[i].index].disso_flag == 1 && particle_list[i].pair_flag != 1 && particle_list[i].update_flag !=1 ){            
            get_parameters_bimol2(species_list, particle_list[i].index, 0, &ka, &D, &sigma, &kd, &gamma_b, &gamma_u);            
            if (ka != 0.0 ){                
                // b = std::min( ALPHA*sigma, b_value + sigma);
                b = ALPHA*sigma;
                // b_value_file << b << std::endl;                
                physparams_bimol2 = make_parameters_bimol2(ka, D, sigma, b, kd, gamma_b, gamma_u);                
                rnd_state = gsl_rng_uniform(gen);
                for ( ii = 0 ; ii < It_coeff_arr.size() ; ii++){
                    if (It_coeff_arr[ii].index == particle_list[i].index){
                        index = ii;
                    }
                }
                It =  It_bimol2(sampled_time, It_coeff_arr[index]);
                Qt =  Qt_bimol2(sampled_time, Qt_coeff_arr[index]);                
                bimol2_state_intermediate = find_state_intermediate(rnd_state, sampled_time, physparams_bimol2, It, Qt);                    
                if (bimol2_state_intermediate == "bound"){                    
                    gaussian_distributing( species_list,  particle_list, i, sampled_time, gen);                    
                }
                else{                    
                    rnd_space = gsl_rng_uniform(gen);//if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
                    sampled_r = draw_space_r_bimol2(rnd_space, sampled_time, physparams_bimol2, space_r_coeff_arr[index], It);
                    // std::cout << sampled_r << std::endl;                    
                    intermediate_unbound(species_list,  particle_list, i, 0, sampled_time, sampled_r, gen);                    
                }
            }
            else{ //ka == 0
                std::cout << __FILE__ << SP << __LINE__ << ": Error: there is a particle with specification beyond the algorithm is developed for." << std::endl;
                std::cout << __FILE__ << SP << __LINE__ << ": Error: this particle should have been treated as a bimolecular transformation." << std::endl;
            }
            
        }
        if (species_list[particle_list[i].index].disso_flag != 1 && particle_list[i].pair_flag != 1 && particle_list[i].update_flag !=1){            
            gaussian_distributing( species_list,  particle_list, i, sampled_time, gen);            
        }
    }

    for (i = 0 ; i < number_of_pairs ; i++){        
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;}         
        get_parameters_bimol1(species_list, particle_list[pair_list[i].p1].index, particle_list[pair_list[i].p2].index, &ka, &D, &sigma);        
        r0 = pair_list[i].dist + species_list[particle_list[pair_list[i].p1].index].radius + species_list[particle_list[pair_list[i].p2].index].radius ;        
        physparams_bimol1 =  make_parameters_bimol1(ka, D, sigma, r0);        
        sampled_r = draw_space_r_bimol1(rnd_space, sampled_time, physparams_bimol1);        
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;}         
        sampled_theta = draw_space_theta_bimol1(rnd_space, sampled_r, sampled_time, physparams_bimol1);        
        bimol1_unbound(species_list, particle_list, pair_list[i].p1, pair_list[i].p2, sampled_time, sampled_r, sampled_theta, gen);        
    }

    for (i = 0 ; i < np ; i++){
        if (particle_list[i].exist_flag == 0){
            particle_list.erase(particle_list.begin()+i);
            np--;
        }
    }
// b_value_file.close();
}

void update_system_monomolecular_reaction(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, std::vector<pair_struct>& pair_list, int pi_monomol,  std::vector<It_coeff_struct> It_coeff_arr, std::vector<Qt_coeff_struct> Qt_coeff_arr, std::vector<space_r_coeff_struct> space_r_coeff_arr, double b_value, double sampled_time, gsl_rng * gen){

    int i = -1;
    int ii = -1;
    int index = -1;
    int np = particle_list.size();
    int number_of_pairs = pair_list.size();
    std::string bimol2_state_intermediate = "--";
    double ka = 0.0;
    double D = 0.0;
    double sigma = 0.0;
    double r0 = 0.0;
    double b = 0.0;
    double kd = 0.0;
    double gamma_b = 0.0;
    double gamma_u = 0.0;
    double rnd_state  = 0.0;
    double rnd_space = 0.0;
    double sampled_r = 0.0;
    double sampled_theta = 0.0;
    double It =  0.0;
    double Qt =  0.0;

    phys_parameters_bimol1 physparams_bimol1;
    phys_parameters_bimol2 physparams_bimol2;
    // std::ofstream  b_value_file;
    // b_value_file.open(DATA_OUTPUT_PATH B_VALUE_FILE_NAME, std::ios_base::app);


    if ( particle_list[pi_monomol].pair_flag != 1){
        monomolecular_reaction( species_list,  particle_list, pi_monomol, sampled_time, gen);
    }

    for (i = 0 ; i < np ; i++ ){
        if (species_list[particle_list[i].index].disso_flag == 1 && particle_list[i].pair_flag < 1){
            get_parameters_bimol2(species_list, particle_list[i].index, 0, &ka, &D, &sigma, &kd, &gamma_b, &gamma_u);
            if (ka != 0.0 ){
                    // b = std::min( ALPHA*sigma, b_value + sigma);
                b = ALPHA*sigma;
                    // b_value_file << b << std::endl;
                    physparams_bimol2 = make_parameters_bimol2(ka, D, sigma, b, kd, gamma_b, gamma_u);
                    rnd_state = gsl_rng_uniform(gen);
                    for ( ii = 0 ; ii < It_coeff_arr.size() ; ii++){
                        if (It_coeff_arr[ii].index == particle_list[i].index){
                            index = ii;
                        }
                    }
                    It =  It_bimol2(sampled_time, It_coeff_arr[index]);
                    Qt =  Qt_bimol2(sampled_time, Qt_coeff_arr[index]);

                    bimol2_state_intermediate = find_state_intermediate(rnd_state, sampled_time, physparams_bimol2, It, Qt);
                if (bimol2_state_intermediate == "bound"){
                    gaussian_distributing( species_list,  particle_list, i, sampled_time, gen);
                }
                else{
                    rnd_space = gsl_rng_uniform(gen);
                    sampled_r = draw_space_r_bimol2(rnd_space, sampled_time, physparams_bimol2, space_r_coeff_arr[index], It);
                     // std::cout << sampled_r << std::endl;
                    intermediate_unbound(species_list, particle_list, i, 0, sampled_time, sampled_r, gen);
                }
            }
            else{ //ka == 0
                std::cout << __FILE__ << SP << __LINE__ << ": Error: there is a particle with specification beyond the algorithm is developed for." << std::endl;
                std::cout << __FILE__ << SP << __LINE__ << ": Error: this particle should have been treated as a bimolecular transformation." << std::endl;
            }
            
        }
        if (species_list[particle_list[i].index].disso_flag != 1 && particle_list[i].pair_flag != 1 && particle_list[i].update_flag != 1){
            gaussian_distributing( species_list,  particle_list, i, sampled_time, gen);
        }
    }
    for (i = 0 ; i < number_of_pairs ; i++){        
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        get_parameters_bimol1(species_list, particle_list[pair_list[i].p1].index, particle_list[pair_list[i].p2].index, &ka, &D, &sigma);
        r0 = pair_list[i].dist + species_list[particle_list[pair_list[i].p1].index].radius + species_list[particle_list[pair_list[i].p2].index].radius ;
        physparams_bimol1 =  make_parameters_bimol1(ka, D, sigma, r0);
        sampled_r = draw_space_r_bimol1(rnd_space, sampled_time, physparams_bimol1);
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        sampled_theta = draw_space_theta_bimol1(rnd_space, sampled_r, sampled_time, physparams_bimol1);
        bimol1_unbound(species_list, particle_list, pair_list[i].p1, pair_list[i].p2, sampled_time, sampled_r, sampled_theta, gen);
        if(pair_list[i].p1 == pi_monomol || pair_list[i].p2 == pi_monomol){
            monomolecular_reaction_bimol_or_monomol_product(species_list,  particle_list, pi_monomol,  gen); // in this case, if the molecule transform to two other, there is a possibility that the particles overlap with the other particle in the pair
        }
    }

    for (i = 0 ; i < np ; i++){
        if (particle_list[i].exist_flag == 0){
            particle_list.erase(particle_list.begin()+i);
            np--;
        }
    }
// b_value_file.close();
}

void update_system_no_reaction(int k , std::vector<Species>& species_list,  std::vector<Particle>& particle_list, std::vector<pair_struct>& pair_list, std::vector<It_coeff_struct> It_coeff_arr, std::vector<Qt_coeff_struct> Qt_coeff_arr, std::vector<space_r_coeff_struct> space_r_coeff_arr, double b_value, double sampled_time, gsl_rng * gen){

    int i = -1;
    int ii = -1;
    int index = -1;
    int np = particle_list.size();
    int number_of_pairs = pair_list.size();
    std::string bimol2_state_intermediate = "--";
    double ka = 0.0;
    double D = 0.0;
    double sigma = 0.0;
    double r0 = 0.0;
    double b = 0.0;
    double kd = 0.0;
    double gamma_b = 0.0;
    double gamma_u = 0.0;
    double rnd_state  = 0.0;
    double rnd_space = 0.0;
    double sampled_r = 0.0;
    double sampled_theta = 0.0;
    double It =  0.0;
    double Qt =  0.0;

    // std::ofstream  b_value_file;
    // b_value_file.open(DATA_OUTPUT_PATH B_VALUE_FILE_NAME, std::ios_base::app);
    phys_parameters_bimol1 physparams_bimol1;
    phys_parameters_bimol2 physparams_bimol2;

    for (i = 0 ; i < np ; i++ ){
        if (species_list[particle_list[i].index].disso_flag == 1 && particle_list[i].pair_flag != 1){
            get_parameters_bimol2(species_list, particle_list[i].index, 0, &ka, &D, &sigma, &kd, &gamma_b, &gamma_u);
            if (ka != 0.0 ){
                    // b = std::min( ALPHA*sigma, b_value + sigma);
                b = ALPHA*sigma;
                    // b_value_file << b << std::endl;
                    physparams_bimol2 = make_parameters_bimol2(ka, D, sigma, b, kd, gamma_b, gamma_u);
                    rnd_state = gsl_rng_uniform(gen);
                    for ( ii = 0 ; ii < It_coeff_arr.size() ; ii++){
                        if (It_coeff_arr[ii].index == particle_list[i].index){
                            index = ii;
                        }
                    }
                    It =  It_bimol2(sampled_time, It_coeff_arr[index]);
                    Qt =  Qt_bimol2(sampled_time, Qt_coeff_arr[index]);
                    bimol2_state_intermediate = find_state_intermediate(rnd_state, sampled_time, physparams_bimol2, It, Qt);
                if (bimol2_state_intermediate == "bound"){
                    gaussian_distributing( species_list,  particle_list, i, sampled_time, gen);
                }
                else{
                    rnd_space = gsl_rng_uniform(gen);
                    sampled_r = draw_space_r_bimol2(rnd_space, sampled_time, physparams_bimol2, space_r_coeff_arr[index], It);
                    intermediate_unbound(species_list,  particle_list, i, 0, sampled_time, sampled_r, gen);
                }
            }
            else{ //ka == 0
                std::cout << __FILE__ << SP << __LINE__ << ": Error: there is a particle with specification beyond the algorithm is developed for." << std::endl;
                std::cout << __FILE__ << SP << __LINE__ << ": Error: this particle should have been treated as a bimolecular transformation." << std::endl;
            }
            
        }
        if (species_list[particle_list[i].index].disso_flag != 1 && particle_list[i].pair_flag != 1 && particle_list[i].update_flag != 1) {
            gaussian_distributing( species_list,  particle_list, i, sampled_time, gen);
        }
    }
    for (i = 0 ; i < number_of_pairs ; i++){
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        get_parameters_bimol1(species_list, particle_list[pair_list[i].p1].index, particle_list[pair_list[i].p2].index, &ka, &D, &sigma);        
        r0 = pair_list[i].dist + species_list[particle_list[pair_list[i].p1].index].radius + species_list[particle_list[pair_list[i].p2].index].radius ;
        physparams_bimol1 =  make_parameters_bimol1(ka, D, sigma, r0);
        sampled_r = draw_space_r_bimol1(rnd_space, sampled_time, physparams_bimol1);
        rnd_space = gsl_rng_uniform(gen); //if (rnd_space < rand_epsilon){rnd_space = 0.0;} 
        sampled_theta = draw_space_theta_bimol1(rnd_space, sampled_r, sampled_time, physparams_bimol1);
        bimol1_unbound(species_list, particle_list, pair_list[i].p1, pair_list[i].p2, sampled_time, sampled_r, sampled_theta, gen);
    }
// b_value_file.close();

}