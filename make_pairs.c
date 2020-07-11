//
//  make_pairs.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/22.
//
//

#include "make_pairs.h"


std::vector<pair_struct> make_pairs(char& three_body_problem, int& number_of_pairs, double& max_diffusion_time_step, double& b_value, std::vector<Species>& species_list, std::vector<Particle>& particle_list, std::vector<particle_particle_diffusion_time_struct>& particle_particle_diffusion_time_array, double temp_max_time_step, double sphere_radius, double H, int k ){

	int pair_index = 0;
	int c = 0;
	int i1 = -1;
	int i2 = -1;
	int np = particle_list.size();
	double Di1 = 0.0;
	double Di2 = 0.0;
	double shell_i1 = 0.0;
	double shell_i2 = 0.0;
	double max_time_step = 0.0;
	double H2 = H*H;
	number_of_pairs = 0;
	three_body_problem = 'F';
	max_diffusion_time_step = 0.0;
	b_value = 0.0;

pair_struct new_element;
std::vector<pair_struct> pair_list;


	if (np > 1){
		while( (pair_index < np*(np-1)/2 )
	 	&&	(particle_particle_diffusion_time_array[pair_index].dist < P_D_F*(species_list[particle_list[particle_particle_diffusion_time_array[pair_index].pi_1].index].radius + species_list[particle_list[particle_particle_diffusion_time_array[pair_index].pi_2].index].radius) ) 
		&& (three_body_problem != 'T' ) ) {
		i1 = particle_particle_diffusion_time_array[pair_index].pi_1;
		i2 = particle_particle_diffusion_time_array[pair_index].pi_2;
		Di1 = species_list[particle_list[i1].index].DiffCoeff;
		Di2 = species_list[particle_list[i2].index].DiffCoeff;
		if( species_list[particle_list[i1].index].disso_flag == 1 || species_list[particle_list[i2].index].disso_flag == 1){
			max_time_step = std::max(particle_particle_diffusion_time_array[pair_index].t_m, particle_list[i2].t_w);
			if (max_time_step < temp_max_time_step){ //comparing the previous time chosen with the total diffuse time to boundary
				// std::cout << __LINE__ << ": the time for I and J gets together or for J to get to boundary is smaller than the total minium time chosen by now, the total minimum time is chosen here then!" <<  std::endl;
				temp_max_time_step = max_time_step;							
			}
			pair_index++;
			// SIGMA = species_list[particle_list[particle_particle_diffusion_time_array[pair_index].pi_1].index].radius + species_list[particle_list[particle_particle_diffusion_time_array[pair_index].pi_2].index].radius;
		}
		else{
		shell_i1 = particle_list[i1].r + species_list[particle_list[i1].index].radius;
		shell_i2 = particle_list[i2].r + species_list[particle_list[i2].index].radius;
		if(particle_list[i1].pair_flag != 1){ // I is not yet in a pair? no
		// std::cout << __LINE__ << ": I is not in a pair!" << std::endl;			
			if(particle_list[i2].pair_flag != 1){ //J is not yet in a pair? no
				// std::cout << __LINE__ << ": J is not in a pair!" << std::endl;
				if(shell_i1 < sphere_radius - INNER_SHELL){ // I is far from boundary? yes
					// std::cout << __LINE__ << ": I is far from boundary!" << std::endl;
					if(shell_i2 < sphere_radius - INNER_SHELL){// J is far from boundary? yes
						// std::cout << __LINE__ << ": J is far from boundary! " << std::endl;
						// std::cout << __LINE__ << ": I and J form a pair!" <<  std::endl;
						// I and J are going to be a new pair
						particle_list[i1].pair_flag = 1;						
						particle_list[i2].pair_flag = 1;
						new_element.p1 = i1;
						new_element.p2 = i2;
						new_element.dist = particle_particle_diffusion_time_array[pair_index].dist;
						pair_list.push_back(new_element);
						number_of_pairs ++;
						pair_index ++;				
						// go back to while loop
					}
					else{// J is close to boundary (but not I)
						// std::cout << __LINE__ << ": J is close to boundary but not I" <<  std::endl;
						//Z: max_time_step = std::min(particle_particle_diffusion_time_array[pair_index].t_m, particle_list[i2].t_w); // comparing the time to the boundary and to diffuse to partner		
						//S:
						max_time_step = std::max(particle_particle_diffusion_time_array[pair_index].t_m, particle_list[i2].t_w);
						if (max_time_step < temp_max_time_step){ //comparing the previous time chosen with the total diffuse time to boundary
							// std::cout << __LINE__ << ": the time for I and J gets together or for J to get to boundary is smaller than the total minium time chosen by now, the total minimum time is chosen here then!" <<  std::endl;
							temp_max_time_step = max_time_step;							
						}
						pair_index ++;
						if (temp_max_time_step < MIN_DIFFUSION_TIME_STEP){
							// std::cout << __LINE__ << ": the time for I and J gets together or for J to get to boundary is smaller than the MIN_DIFFUSION_TIME_STEP, therefore we adjust all particles!" <<  std::endl;
							three_body_problem = 'T'; //Case6
						}
					}
				}
				else{ // I is close to boundary (we don't know about J)
					// std::cout << __LINE__ << ": I is close to boundary but we don't know about J!"  << std::endl;
					//Z: max_time_step = std::min(std::min(particle_particle_diffusion_time_array[pair_index].t_m, particle_list[i1].t_w), particle_list[i2].t_w);
					//S:					
					max_time_step = std::max(particle_particle_diffusion_time_array[pair_index].t_m, particle_list[i1].t_w);
					if (max_time_step < temp_max_time_step){
					// std::cout << __LINE__ << ": the time for I and J gets together or for I to get to boundary is smaller than the total minium time chosen by now, the total minimum time is chosen here then!" <<  std::endl;			
						temp_max_time_step = max_time_step;
					}
					pair_index ++;
					if (temp_max_time_step < MIN_DIFFUSION_TIME_STEP){
						// std::cout << __LINE__ << ": the time for I and J gets together or for J to get to boundary is smaller than the MIN_DIFFUSION_TIME_STEP, therefore we adjust all particles!" <<  std::endl;						
						three_body_problem = 'T'; //Case 7
					}
				}
			}
			else{//J is in a pair but not I
				// std::cout << __LINE__ << ": J is in a pair!" << std::endl;
				if (shell_i1 > sphere_radius - INNER_SHELL){ // is I close to the boundary? yes
					// std::cout << __LINE__ << ": I close to boundary!" << std::endl;
					//Z: max_time_step = std::min(particle_particle_diffusion_time_array[pair_index].t_m, particle_list[i1].t_w);
					//S:
					max_time_step = particle_particle_diffusion_time_array[pair_index].t_m;
					if (max_time_step < temp_max_time_step){
						// std::cout << __LINE__ << ": the time for I to get to J or for I to get to boundary is smaller than the total minimum time chosen by now, the total minimum time is chosen here!" << std::endl;			
						temp_max_time_step = max_time_step;
					}
					pair_index ++;
					if (temp_max_time_step < MIN_DIFFUSION_TIME_STEP){
						// std::cout << __LINE__ << ": the time for I to get to J or for I to get to boundary is smaller than MIN_DIFFUSION_TIME_STEP, therefore we adjust all particles!" << std::endl;
						three_body_problem = 'T'; //Case3
					}
				}
				else{ // I is not close to the boundary
					// std::cout << __LINE__ << ": I is not close to boundary!" << std::endl;
					//Z & S:
					max_time_step = particle_particle_diffusion_time_array[pair_index].t_m;				
					if (max_time_step < temp_max_time_step){
						temp_max_time_step = max_time_step;					
					}
					pair_index ++;
					if (temp_max_time_step < MIN_DIFFUSION_TIME_STEP){
						three_body_problem = 'T'; //Case 5
					}

				}
			}
		}
		else{ // I is in a pair but we don't know about J
			// std::cout << __LINE__ << ": I is in a pair but we don't know about J" << std::endl;
			if (particle_list[i2].pair_flag == 1){
				// std::cout << __LINE__ << ": I and J are in a pair already!" << std::endl;
				max_time_step = particle_particle_diffusion_time_array[pair_index].t_m;				
				if (max_time_step < temp_max_time_step){
					temp_max_time_step = max_time_step;					
				}
				pair_index ++;
				if (temp_max_time_step < MIN_DIFFUSION_TIME_STEP){
					three_body_problem = 'T'; //Case 1
				}
			}
			else{ //J is not in a pair!
				// std::cout << __LINE__ << ": J is not in a pair!" << std::endl;
				if (shell_i2 > sphere_radius - INNER_SHELL){ // J is close to boundary!
					// std::cout << __LINE__ << ": J is close to boundary!"<< std::endl;

					//Z: max_time_step = std::min(particle_particle_diffusion_time_array[pair_index].t_m, particle_list[i2].t_w);
					//S:
					max_time_step = particle_particle_diffusion_time_array[pair_index].t_m;
					if (max_time_step < temp_max_time_step){
						// std::cout << __LINE__ << ": the time for J to get to I or get to boundary is smaller than the total minimum time chosen by now, the new total minimum time is chosen here!" << std::endl;
						temp_max_time_step = max_time_step;
					}
					pair_index ++;
					if (temp_max_time_step < MIN_DIFFUSION_TIME_STEP){
						// std::cout << __LINE__ << ": the time for J to get to I or get to boundary is smaller than MIN_DIFFUSION_TIME_STEP, we need to adjust all particles!" << std::endl;
						three_body_problem = 'T'; //Case2
					}
				}
				else{ // J is not close to the boundary
					// std::cout << __LINE__ << ": J is not close to the boundary!" <<std::endl;
					//Z & S:
					max_time_step = particle_particle_diffusion_time_array[pair_index].t_m;
					if (max_time_step < temp_max_time_step){
						// std::cout << __LINE__ <<": the time for J to get to I is smaller than the total minimum time chosen by now, the new total minimum time is chosen here!" << std::endl;
						temp_max_time_step = max_time_step;					
					}
					pair_index ++;
					if (temp_max_time_step < MIN_DIFFUSION_TIME_STEP){
						// std::cout << __LINE__ << ": the time for J to get to I is smaller than MIN_DIFFUSION_TIME_STEP, we need to adjust all particles!" << std::endl;
						three_body_problem = 'T'; //Case4
					}
				}
			}

		}
	}	
	}

	}
	


	max_diffusion_time_step = temp_max_time_step;
	
    return pair_list;	

}

