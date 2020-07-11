//
//  bimolecular_reaction_type_two.c
//  
//
//  Created by Zahedeh Bashardanesh on 09/05/14.
//
//
 

#include "bimolecular_reaction_type_two.h"

void intermediate_unbound(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_bimol2, int J_disso, double sampled_time, double sampled_r, gsl_rng * gen){

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

	struct Particle new_particle;

	DC = species_list[particle_list[pi_bimol2].index].DiffCoeff;
	DA = species_list[sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product1.c_str(), species_list)].DiffCoeff ;
	DB = species_list[sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product2.c_str(), species_list)].DiffCoeff;
	

	DA_p_DB = DA + DB;
	DA_p_DB_r = 1.0/DA_p_DB;

	alpha = sqrt((DB/DA)*DC*DA_p_DB_r);
	beta = sqrt((DA/DB)*DC*DA_p_DB_r);
	
	// alpha = DB*DA_p_DB_r;
	// beta = DA*DA_p_DB_r;

	alpha_p_beta = alpha + beta;
	alpha_p_beta_r = 1.0/alpha_p_beta;


	dx = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));                
    dy = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));
    dz = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));


    //Initial position for the center of diffusion vector follows the R = alpha*rA + beta*rB, with considering rA and rB are equal to rC
	x_R =  (alpha_p_beta) * particle_list[pi_bimol2].x + dx;
	y_R =  (alpha_p_beta) * particle_list[pi_bimol2].y + dy;
	z_R =  (alpha_p_beta) * particle_list[pi_bimol2].z + dz;

    x_r = gsl_ran_gaussian(gen, 1);
    y_r = gsl_ran_gaussian(gen, 1);
    z_r = gsl_ran_gaussian(gen, 1);
    r_r = x_r * x_r + y_r * y_r + z_r * z_r;
    r_r = sqrt(r_r); 
    r_r_r = 1.0/r_r;
    x_r = sampled_r * x_r * r_r_r;
    y_r = sampled_r * y_r * r_r_r;
    z_r = sampled_r * z_r * r_r_r;

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
	new_particle.phi = atan(new_particle.y/new_particle.x);
	new_particle.theta = acos(new_particle.z/new_particle.r);	
	new_particle.index = sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product2.c_str(), species_list);
	new_particle.update_flag++;
	particle_list.push_back(new_particle);
	

	particle_list[pi_bimol2].x = x_A;
	particle_list[pi_bimol2].y = y_A;
	particle_list[pi_bimol2].z = z_A;
	particle_list[pi_bimol2].r = sqrt(particle_list[pi_bimol2].x*particle_list[pi_bimol2].x + particle_list[pi_bimol2].y*particle_list[pi_bimol2].y + particle_list[pi_bimol2].z*particle_list[pi_bimol2].z);
	if(particle_list[pi_bimol2].x < 0.0){
		particle_list[pi_bimol2].phi = atan(particle_list[pi_bimol2].y/particle_list[pi_bimol2].x) + Pi;	
	}
	else{
		particle_list[pi_bimol2].phi = atan(particle_list[pi_bimol2].y/particle_list[pi_bimol2].x);
	}
	particle_list[pi_bimol2].theta = acos(particle_list[pi_bimol2].z/particle_list[pi_bimol2].r);
	particle_list[pi_bimol2].index = sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product1.c_str(), species_list);
	particle_list[pi_bimol2].update_flag ++;


}



void separation(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_bimol2, int J_disso, double sampled_time, double sep_dist, gsl_rng * gen){

// std::cout << "2" << std::endl;
	int last_element = particle_list.size();



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

	double b = sep_dist;

	struct Particle new_particle;


	DC = species_list[particle_list[pi_bimol2].index].DiffCoeff;
	DA = species_list[sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product1.c_str(), species_list)].DiffCoeff;
	DB = species_list[sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product2.c_str(), species_list)].DiffCoeff;
	


	DA_p_DB = DA + DB;
	DA_p_DB_r = 1.0/DA_p_DB;

	alpha = sqrt((DB/DA)*DC*DA_p_DB_r);
	beta = sqrt((DA/DB)*DC*DA_p_DB_r);
	
	// alpha = DB*DA_p_DB_r;
	// beta = DA*DA_p_DB_r;

	alpha_p_beta = alpha + beta;
	alpha_p_beta_r = 1.0/alpha_p_beta;
	

	dx = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));                
    dy = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));
    dz = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));



    //The locaiton of particle C considered to be the initial location for the center of diffusion vector
    // x_R =  particle_list[pi_bimol2].x  + dx;
    // y_R =  particle_list[pi_bimol2].y  + dy;
    // z_R =  particle_list[pi_bimol2].z  + dz;


    //Initial position for the center of diffusion vector follows the R = alpha*rA + beta*rB, with considering rA and rB are equal to rC
	x_R =  (alpha_p_beta)*particle_list[pi_bimol2].x + dx;
	y_R =  (alpha_p_beta)*particle_list[pi_bimol2].y + dy;
	z_R =  (alpha_p_beta)*particle_list[pi_bimol2].z + dz;

    x_r = gsl_ran_gaussian(gen, 1);
    y_r = gsl_ran_gaussian(gen, 1);
    z_r = gsl_ran_gaussian(gen, 1);
    r_r = x_r * x_r + y_r * y_r + z_r * z_r;
    r_r = sqrt(r_r);
    r_r_r = 1.0/r_r;

    x_r = b * x_r * r_r_r;
    y_r = b * y_r * r_r_r;
    z_r = b * z_r * r_r_r;


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
	new_particle.index = sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product2.c_str(), species_list);
	new_particle.update_flag++;
	particle_list.push_back(new_particle);

	particle_list[pi_bimol2].x = x_A;
	particle_list[pi_bimol2].y = y_A;
	particle_list[pi_bimol2].z = z_A;
	particle_list[pi_bimol2].r = sqrt(particle_list[pi_bimol2].x*particle_list[pi_bimol2].x + particle_list[pi_bimol2].y*particle_list[pi_bimol2].y + particle_list[pi_bimol2].z*particle_list[pi_bimol2].z);
		if(particle_list[pi_bimol2].x < 0.0){
		particle_list[pi_bimol2].phi = atan(particle_list[pi_bimol2].y/particle_list[pi_bimol2].x) + Pi;	
	}
	else{
		particle_list[pi_bimol2].phi = atan(particle_list[pi_bimol2].y/particle_list[pi_bimol2].x);
	}
	particle_list[pi_bimol2].theta = acos(particle_list[pi_bimol2].z/particle_list[pi_bimol2].r);
	particle_list[pi_bimol2].index = sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product1.c_str(), species_list);
	particle_list[pi_bimol2].update_flag ++;

}

void exit_bound(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int pi_bimol2, double sampled_time, gsl_rng * gen){

// the particle that has been chosen for dissociation event is going to transform while it's been in the bound state
// at this state it can transform to one other particle through different channels (monomolecular channels) or to other
// particles through different channels (bimolecular channels) or to decay, we have to see which channel is chosen 

	gaussian_distributing(species_list,  particle_list, pi_bimol2, sampled_time, gen);
	monomolecular_reaction_bimol_or_monomol_product(species_list, particle_list, pi_bimol2, gen);

}



void exit_unbound(std::vector<Species>& species_list, std::vector<Particle>& particle_list, int pi_bimol2, int J_disso, double sampled_time, double sampled_r, gsl_rng * gen){

// the possiblity of pi_bimol2 can go to more than one (rev)dissociation channel is not take into account therefore J_disso = 0;
	J_disso = 0;
	int last_element = particle_list.size();
	int i = 0;
	int transf_partner_1 = -1;
	int transf_partner_2 = -1;
	int partner_index = -1;
	int switch_01 = 0;
	int A_index = -1;
	int B_index = -1;

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

	
	double a = 0.0;
	double a_sum = 0.0;
	double a_sum_r = 0.0;
	double b = 0.0;
	double b_sum = 0.0;
	double b_sum_r = 0.0;
	double ab_sum = 0.0;
	double ab_sum_r = 0.0;
	double u = 0.0;

	struct Particle new_particle;



	DC = species_list[particle_list[pi_bimol2].index].DiffCoeff;
	DA = species_list[sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product1.c_str(), species_list)].DiffCoeff ;
	DB = species_list[sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product2.c_str(), species_list)].DiffCoeff;
	

	
	DA_p_DB = DA + DB;
	DA_p_DB_r = 1.0/DA_p_DB;

	alpha = sqrt((DB/DA)*DC*DA_p_DB_r);
	beta = sqrt((DA/DB)*DC*DA_p_DB_r);
	

	// alpha = DB*DA_p_DB_r;
	// beta = DA*DA_p_DB_r;
	
	alpha_p_beta = alpha + beta;
	alpha_p_beta_r = 1.0/alpha_p_beta;


	dx = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));                
    dy = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));
    dz = gsl_ran_gaussian(gen, sqrt(2*DC*sampled_time));

    x_R =  (alpha_p_beta) * particle_list[pi_bimol2].x + dx;
	y_R =  (alpha_p_beta) * particle_list[pi_bimol2].y + dy;
	z_R =  (alpha_p_beta) * particle_list[pi_bimol2].z + dz;

    x_r = gsl_ran_gaussian(gen, 1);
    y_r = gsl_ran_gaussian(gen, 1);
    z_r = gsl_ran_gaussian(gen, 1);
    r_r = x_r * x_r + y_r * y_r + z_r * z_r;
    r_r = sqrt(r_r);
    r_r_r = 1.0/r_r;
    x_r = sampled_r * x_r * r_r_r;
    y_r = sampled_r * y_r * r_r_r;
    z_r = sampled_r * z_r * r_r_r;

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
	new_particle.phi = atan(new_particle.y/new_particle.x);
	new_particle.theta = acos(new_particle.z/new_particle.r);
	new_particle.index = sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product2.c_str(), species_list);
	new_particle.update_flag++;
	particle_list.push_back(new_particle);
	

	particle_list[pi_bimol2].x = x_A;
	particle_list[pi_bimol2].y = y_A;
	particle_list[pi_bimol2].z = z_A;
	particle_list[pi_bimol2].r = sqrt(particle_list[pi_bimol2].x*particle_list[pi_bimol2].x + particle_list[pi_bimol2].y*particle_list[pi_bimol2].y + particle_list[pi_bimol2].z*particle_list[pi_bimol2].z);
	particle_list[pi_bimol2].phi = atan(particle_list[pi_bimol2].y/particle_list[pi_bimol2].x);
	particle_list[pi_bimol2].theta = acos(particle_list[pi_bimol2].z/particle_list[pi_bimol2].r);
	particle_list[pi_bimol2].update_flag++;
	particle_list[pi_bimol2].index = sttoint(species_list[particle_list[pi_bimol2].index].disso[J_disso].product1.c_str(), species_list);
	

	A_index = particle_list[pi_bimol2].index;
	B_index = particle_list[last_element].index;

	a_sum = species_list[A_index].total_transform_rate;
	b_sum = species_list[B_index].total_transform_rate;

	ab_sum = a_sum + b_sum;
	ab_sum_r = 1.0/ab_sum;

	u = gsl_rng_uniform(gen);
	if (u < a_sum*ab_sum_r){ // particle A is going to transform
		monomolecular_reaction_bimol_or_monomol_product(species_list, particle_list, pi_bimol2, gen);
	}
	else{ //particle B is going to transform
		monomolecular_reaction_bimol_or_monomol_product(species_list, particle_list, last_element, gen);
	}
}

