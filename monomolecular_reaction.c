//
//  monomolecular_reaction.c
//  
//
//  Created by Zahedeh Bashardanesh on 09/23/14.
//
//

 
#include "monomolecular_reaction.h"

void monomolecular_reaction(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_monomol, double sampled_time, gsl_rng * gen){



	gaussian_distributing(species_list,  particle_list, pi_monomol, sampled_time, gen);
	monomolecular_reaction_bimol_or_monomol_product(species_list, particle_list, pi_monomol, gen);
	
}


void monomolecular_reaction_bimol_or_monomol_product(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_monomol,  gsl_rng * gen){

	int pi_monomol_index = particle_list[pi_monomol].index;
	int i = 0;
	int j = 0;
	int transfer_to_index = -1;
	int J = -1;
	int monomol_flag = -1;
	int bimol_flag = -1;
	double u = 0.0;
	double a_sum_r = 0.0;
	double a = 0.0;
	
		if (species_list[pi_monomol_index].trans_monomol_flag != 0) {
			if (species_list[pi_monomol_index].trans_bimol_flag == 1){ // both mono- and bi-molecular transformation can happen
				u = gsl_rng_uniform(gen);
				a_sum_r = 1.0/species_list[pi_monomol_index].total_transform_rate;
				for (i = 0 ; i < species_list[pi_monomol_index].trans_monomol_channel ; i++){
					a = a + species_list[pi_monomol_index].trans_monomol[i].trans_monomol_rate;
					if (u < a*a_sum_r){
						transfer_to_index = i;
						i = species_list[pi_monomol_index].trans_monomol_channel;
						J = species_list[pi_monomol_index].trans_bimol_channel;
						monomol_flag = 1;
					}
				}
				for (j = J ; j < species_list[pi_monomol_index].trans_bimol_channel ; j ++){
					a = a + species_list[pi_monomol_index].trans_bimol[j].trans_bimol_rate;
					if (u < a*a_sum_r){
						transfer_to_index = j;
						j = species_list[pi_monomol_index].trans_bimol_channel;
						bimol_flag = 1;
					}
				}
				if (monomol_flag == 1){ // monomolecular transformation happens
					if ( species_list[pi_monomol_index].trans_monomol[transfer_to_index].product == "0"){//decay
						particle_list[pi_monomol].exist_flag = 0;
					}
					else{ //transform
						particle_list[pi_monomol].index = sttoint(species_list[pi_monomol_index].trans_monomol[transfer_to_index].product.c_str(), species_list);
						particle_list[pi_monomol].update_flag ++;
					}
				}
				else if (bimol_flag == 1){ //bimolecular transformation happens
					monomolecular_reaction_bimol_product(species_list, particle_list, pi_monomol, transfer_to_index, gen);
				}
				else{
					std::cout << __LINE__  << ": error: no transformation happening for this species!" << std::endl;
				}

			}
			else{//just monomolecular transformation
			u = gsl_rng_uniform(gen);
			a = species_list[pi_monomol_index].trans_monomol[0].trans_monomol_rate;
			i++;
			a_sum_r = 1.0/species_list[pi_monomol_index].total_transform_rate;
			while ( u > a*a_sum_r){
				a = a + species_list[pi_monomol_index].trans_monomol[i].trans_monomol_rate;
				i++;
			}
			if ( species_list[pi_monomol_index].trans_monomol[i-1].product == "0"){//decay
				particle_list[pi_monomol].exist_flag = 0;
			}
			else{ //transform
				particle_list[pi_monomol].index = sttoint(species_list[pi_monomol_index].trans_monomol[i-1].product.c_str(), species_list);
				particle_list[pi_monomol].update_flag ++;
			}

		}
	}
	else if(species_list[pi_monomol_index].trans_bimol_flag != 0){//just bimolecular transformation
		u = gsl_rng_uniform(gen);
		a = species_list[pi_monomol_index].trans_bimol[0].trans_bimol_rate;
		
		a_sum_r = 1.0/species_list[pi_monomol_index].total_transform_rate;
		while ( u > a*a_sum_r){
			i++;
			a = a + species_list[pi_monomol_index].trans_bimol[i].trans_bimol_rate;
		}
		monomolecular_reaction_bimol_product(species_list, particle_list, pi_monomol, i, gen);
	}

	else{//no transformation 
		std::cout << __LINE__ <<  ": error: no transformation happening for this species!" << std::endl;
	}
}



void monomolecular_reaction_bimol_product(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_monomol,  int transfer_to_index, gsl_rng * gen){


	int last_element = particle_list.size();
	int J_disso = 0; // there is only the possiblity that one partcile reversibly dissociate to one pair
	double DA = 0.0;
	double DB = 0.0;
	double DA_p_DB = 0.0;
	double DA_p_DB_r = 0.0;
	double DC = 0.0;
	double alpha = 0.0;
	double beta = 0.0;
	double alpha_p_beta = 0.0;
	double alpha_p_beta_r = 0.0;
	double dx = 0.0;
	double dy  = 0.0;
	double dz = 0.0;
	double x_R = 0.0;
	double y_R = 0.0;
	double z_R = 0.0;

	double x_A = 0.0;
	double y_A = 0.0;
	double z_A = 0.0;	
	double x_B = 0.0;
	double y_B = 0.0;
	double z_B = 0.0;	

	double x_r = 0.0;
	double y_r = 0.0;
	double z_r = 0.0;
	double r_r = 0.0;	
	double r_r_r = 0.0;
	double radii = 0.0;
	double radius_A = 0.0;
	double radius_B = 0.0;
	struct Particle new_particle;

	DC = species_list[particle_list[pi_monomol].index].DiffCoeff;
	DA = species_list[sttoint(species_list[particle_list[pi_monomol].index].trans_bimol[transfer_to_index].product1.c_str(), species_list)].DiffCoeff ;
	DB = species_list[sttoint(species_list[particle_list[pi_monomol].index].trans_bimol[transfer_to_index].product2.c_str(), species_list)].DiffCoeff;
	radius_A = species_list[sttoint(species_list[particle_list[pi_monomol].index].trans_bimol[transfer_to_index].product1.c_str(), species_list)].radius;
	radius_B = species_list[sttoint(species_list[particle_list[pi_monomol].index].trans_bimol[transfer_to_index].product2.c_str(), species_list)].radius;
	radii = radius_A + radius_B;

	DA_p_DB = DA + DB;
	DA_p_DB_r = 1.0/DA_p_DB;

	// alpha = sqrt((DB/DA)*DC*DA_p_DB_r);
	// beta = sqrt((DA/DB)*DC*DA_p_DB_r);
	
	alpha = DB*DA_p_DB_r;
	beta = DA*DA_p_DB_r;

	// alpha_p_beta = sqrt(alpha) + sqrt(beta);
	alpha_p_beta = alpha + beta;
	alpha_p_beta_r = 1.0/alpha_p_beta;


    
	x_R =  particle_list[pi_monomol].x ;
	y_R =  particle_list[pi_monomol].y ;
	z_R =  particle_list[pi_monomol].z ;
	

    x_r = gsl_ran_gaussian(gen, 1);
    y_r = gsl_ran_gaussian(gen, 1);
    z_r = gsl_ran_gaussian(gen, 1);
    r_r = x_r * x_r + y_r * y_r + z_r * z_r;
    r_r = sqrt(r_r); 
    r_r_r = 1.0/r_r;
    x_r = radii* x_r * r_r_r;
    y_r = radii* y_r * r_r_r;
    z_r = radii* z_r * r_r_r;

	x_A = (x_R - beta * x_r)*alpha_p_beta_r;
	y_A = (y_R - beta * y_r)*alpha_p_beta_r;
	z_A = (z_R - beta * z_r)*alpha_p_beta_r;

	x_B = (x_R + alpha * x_r)*alpha_p_beta_r;
	y_B = (y_R + alpha * y_r)*alpha_p_beta_r;
	z_B = (z_R + alpha * z_r)*alpha_p_beta_r;

	new_particle.x = x_B;
	new_particle.y = y_B;
	new_particle.z = z_B;
	new_particle.r = sqrt(new_particle.x*new_particle.x + new_particle.y*new_particle.y + new_particle.z*new_particle.z);
	if(new_particle.x < 0.0){
		new_particle.phi = atan(new_particle.y/new_particle.x) + Pi;	
	}
	else{
		new_particle.phi = atan(new_particle.y/new_particle.x);
	}
	new_particle.theta = acos(new_particle.z/new_particle.r);	
	new_particle.index = sttoint(species_list[particle_list[pi_monomol].index].trans_bimol[transfer_to_index].product2.c_str(), species_list);
	particle_list.push_back(new_particle);
	particle_list[last_element].update_flag ++;


	particle_list[pi_monomol].x = x_A;
	particle_list[pi_monomol].y = y_A;
	particle_list[pi_monomol].z = z_A;
	particle_list[pi_monomol].r = sqrt(particle_list[pi_monomol].x*particle_list[pi_monomol].x + particle_list[pi_monomol].y*particle_list[pi_monomol].y + particle_list[pi_monomol].z*particle_list[pi_monomol].z);
	if(particle_list[pi_monomol].x < 0.0){
		particle_list[pi_monomol].phi = atan(particle_list[pi_monomol].y/particle_list[pi_monomol].x) + Pi;	
	}
	else{
		particle_list[pi_monomol].phi = atan(particle_list[pi_monomol].y/particle_list[pi_monomol].x);
	}
	particle_list[pi_monomol].theta = acos(particle_list[pi_monomol].z/particle_list[pi_monomol].r);
	particle_list[pi_monomol].index = sttoint(species_list[particle_list[pi_monomol].index].trans_bimol[transfer_to_index].product1.c_str(), species_list);
	particle_list[pi_monomol].update_flag ++;

}






