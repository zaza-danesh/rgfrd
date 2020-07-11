//
//  main_rgfrd.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/21.
//
//
 
// #ifdef __linux__
//     #include  LIB_PATH 
// #else 
//     #include LIB_PATH 
// #endif

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif



#include "mystructure.h"
#include "get_system.h"
#include "get_species.h"
#include "get_reactions.h"
#include "calculate_distances.h"
#include "initial_distribution.h"
#include "adjust_particles.h"
#include "make_pairs.h"
#include "particle_particle_diffusion.h"
#include "particle_boundary_diffusion.h"
#include "check_max_time_step.h"
#include "draw_event_time.h"
#include "update_system.h"
#include "FindRoots.h"
#include "It_coeff.h"
#include "Qt_coeff.h"
#include "Jt_coeff.h"
#include "space_r_coeff.h"
#include "reflection.h"
#include "save_data.h"


int main(int argc, char *argv[]){


	int i = -1;
	int k = -1;
	int np = 0;
	int ns = -1;
	int pair_index = -1;
	int number_of_pairs = -1;
	int short_time_counter = 0;
	int event_flag = -1;
	double t = 0.0;
	double final_time = atof(argv[2]);
	double sampled_time = 0.0;
	double b_value = 0.0;
	double max_diffusion_time_step = 0.0;
    double seed = atof(argv[1]);
	double H = 2.0;
	double R_curv = 5;
	char negative_distances = 'T';
	char three_body_problem = 'F';
	double ref_flag = 1;
	std::string event = "--";

	gsl_rng * gen;  //generator
    const gsl_rng_type * T;
    T = gsl_rng_default; //default is mt19937
    gen = gsl_rng_alloc (T);
    gsl_rng_set(gen, seed);


    double S_init = 10;
	double P_init = 10;
	double E1_init = 1;
	double E2_init = 1;
	double SE1_init = 0;
	double PE2_init = 0;
	double SP_tot = S_init + P_init;

	double delta_t = 1e-2;
	double sigma = 1e-2;
	double Rp = 2*sigma; //2sigma
    double diffusion_coefficient =  sigma*sigma/delta_t;
    double Dp = 2*diffusion_coefficient; //2D
    double kD = 4*Pi*Rp*Dp;
    double sphere_radius = .1;
    double Vol = 4.0*Pi*(sphere_radius*sphere_radius*sphere_radius)/3.0;
    double ka= (atof(argv[3]) - 1.0)*kD;
    double kd = sigma*sigma/diffusion_coefficient;
    double k1 = sigma*sigma/diffusion_coefficient;
    double k2 = sigma*sigma/diffusion_coefficient;


        
	double Number_of_particles = S_init + P_init + E1_init + E2_init + SE1_init + PE2_init;
	np = Number_of_particles;
	std::ofstream species_file;
	species_file.open(DATA_OUTPUT_PATH "species.txt");
	species_file << "Number_of_species: 6 " << std::endl; 
	species_file << "Name  radius DiffCoeff initial_quantity " << std::endl;
	species_file << "S" 	<<  SP << sigma << SP << diffusion_coefficient << SP << S_init << std::endl;
	species_file << "P" 	<<  SP << sigma << SP << diffusion_coefficient << SP << P_init << std::endl;
	species_file << "E1" 	<<  SP << sigma << SP << diffusion_coefficient << SP << E1_init << std::endl;
	species_file << "E2" 	<<  SP << sigma << SP << diffusion_coefficient << SP << E2_init << std::endl;
	species_file << "SE1" 	<<  SP << sigma << SP << diffusion_coefficient << SP << SE1_init << std::endl;
	species_file << "PE2" 	<<  SP << sigma << SP << diffusion_coefficient << SP << PE2_init << std::endl;

	species_file.close();

	std::vector < std::vector <double> > dij_matrix (np, std::vector<double>(np));
	std::vector<Particle> particle_list;

	std::vector<particle_particle_diffusion_time_struct> particle_particle_diffusion_time_array;
	std::vector<particle_boundary_diffusion_time_struct> particle_boundary_diffusion_time_array;
	std::vector<pair_struct> pair_list;
	event_struct reaction_event;


	std::vector<Species> species_list;
	std::vector<Roots> species_roots;
	std::vector<It_coeff_struct> It_coeff_arr;
	std::vector<Qt_coeff_struct> Qt_coeff_arr;
	std::vector<Jt_coeff_struct> Jt_coeff_arr;
	std::vector<space_r_coeff_struct> space_r_coeff_arr;	

	get_reactions(get_number_of_species(), species_list, ka, kd, k1, ka, kd, k2);
	particle_list = initial_distribution(species_list, dij_matrix, np, sphere_radius, gen, k);

   	species_roots = Find_Roots(species_list, b_value);
   	It_coeff_arr =	It_coeff(species_list, species_roots, b_value);
   	Qt_coeff_arr =	Qt_coeff(species_list, species_roots, b_value);
   	Jt_coeff_arr =	Jt_coeff(species_list, species_roots, b_value);
   	space_r_coeff_arr =	space_r_coeff(species_list, species_roots, b_value);
	
	double t1 = clock();
	// double kk = 1;
	// int SE1 = 0;
	// int PE2 = 0;
    while( t < final_time){

    	// if (t > kk  *atof(argv[4]) ){
    	// 	// std::cout << kk << std::endl;
    	// 	std::cout << t << std::endl;
    	// 	kk++;
    	// }

		k++;
		sampled_time = 0.0;
   		number_of_pairs = 0;
   		max_diffusion_time_step = 0.0; 
    	if ( k != 0){
    		if (np != particle_list.size()){
    			np = particle_list.size();
    			dij_matrix.resize(np);
    			for (i = 0 ; i < dij_matrix.size() ; i++){
    				dij_matrix[i].resize(np);	
    			}    			
    		}		
    		if (np > 1){
    			negative_distances = calculate_distances(species_list, particle_list, dij_matrix, np, k);    			
    		}
    		for ( i = 0 ; i < np ; i++)	{
    			particle_list[i].update_flag = 0;
    			particle_list[i].pair_flag = 0;
    			particle_list[i].exist_flag = 1;
    			particle_list[i].reflect_flag = 0;
    		}
    	}
    	else{
    		negative_distances = 'F';
    	}    	
    	while (negative_distances == 'T'){
    		adjust_particles(species_list, particle_list, sphere_radius, k );
    		negative_distances = calculate_distances(species_list, particle_list, dij_matrix, np, k);   		
    	}
    	if (np>1){
    		particle_particle_diffusion_time_array = particle_particle_diffusion(species_list, particle_list, dij_matrix, H, k );	
    	}
    	particle_boundary_diffusion_time_array = particle_boundary_diffusion(species_list, particle_list, dij_matrix, H, sphere_radius, R_curv, k );
    	if (np>1){
    		pair_list = make_pairs(three_body_problem, number_of_pairs, max_diffusion_time_step, b_value, species_list, particle_list,  particle_particle_diffusion_time_array, particle_boundary_diffusion_time_array[0].t_w, sphere_radius, H, k );
    	}    	
    	if (three_body_problem == 'T'){
    		adjust_particles(species_list, particle_list, sphere_radius, k );
    	}
    	else{
    		if (check_max_time_step(max_diffusion_time_step, short_time_counter) != 1 ){
    			reaction_event = draw_event_time(species_list, particle_list, pair_list, number_of_pairs,  It_coeff_arr, Qt_coeff_arr, b_value, gen, k);
				if ( reaction_event.t_r < max_diffusion_time_step ){
					event = reaction_event.event;
					sampled_time = reaction_event.t_r;
					event_flag = 1;
				}
				else{
					event = "no_reaction";
					sampled_time = max_diffusion_time_step;
					event_flag = 0;
				}

				if ( event_flag == 1 ) {
					if ( event == "bimol1" ){
						t = t + sampled_time;
						update_system_bimolecular_reaction_type_one(species_list, particle_list, pair_list, reaction_event.pair_index,  It_coeff_arr, Qt_coeff_arr, space_r_coeff_arr,  b_value, sampled_time, gen);
						reflection(species_list, particle_list, sphere_radius, ref_flag, final_time, k);						
					}//else if bimol1
					else if ( event == "bimol2" ){
						t = t + sampled_time;
						update_system_bimolecular_reaction_type_two(k, species_list, particle_list, pair_list, reaction_event.pi_bimol2,  It_coeff_arr, Qt_coeff_arr, Jt_coeff_arr, space_r_coeff_arr,  b_value, sampled_time, gen);
						reflection(species_list, particle_list, sphere_radius, ref_flag,  final_time, k);
					}//elsif bimol1
					else if ( event == "monomol" ){
						t = t + sampled_time;
						update_system_monomolecular_reaction(species_list, particle_list, pair_list, reaction_event.pi_monomol,  It_coeff_arr, Qt_coeff_arr, space_r_coeff_arr,  b_value, sampled_time, gen);
						reflection(species_list, particle_list, sphere_radius, ref_flag,  final_time, k);				
					} //elsif monomol
					else {
						std::cout << __FILE__ << SP << __LINE__ <<  ": Error: no event chosen." << std::endl;	
					}
				}
				else {
					t = t + sampled_time; 
					update_system_no_reaction(k, species_list, particle_list, pair_list,  It_coeff_arr, Qt_coeff_arr, space_r_coeff_arr,  b_value, sampled_time, gen);
					reflection(species_list, particle_list, sphere_radius, ref_flag,  final_time, k);

				}//else, no reaction happens
			}//if check_max_time
    		else{
    			t = t + final_time;
    		}
    		// SE1 = 0;
    		// PE2 = 0;
    		// for (i = 0 ; i < particle_list.size() ; i++){
    		// 	if (particle_list[i].index == 5){
    		// 		SE1++;
    		// 		// std::cout << "SE1" << std::endl;
    		// 	}
    		// 	else if (particle_list[i].index == 6){
    		// 		PE2++;
    		// 		// std::cout << "PE2" << std::endl;
    		// 	}
    		// }
    		// if (SE1 != 0 || PE2 != 0){
    		// 	std::cout << SE1 << SP << PE2 << std::endl;
    		// }

    	}
    }//while
    double t2 = clock();

    std::cout << (t2 - t1)/CLOCKS_PER_SEC << std::endl;
}//main
