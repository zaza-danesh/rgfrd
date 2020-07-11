//
//  bimolecular_reaction_type_one.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/07/14.
//
//
 

#include "bimolecular_reaction_type_one.h"

void bimol1_bound(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int p1, int p2, double sampled_time, gsl_rng * gen){


	int index = 0;
	int p1_index = -1;
	int p2_index = -1;

	double D_A = 0.0;
	double D_B = 0.0;
	double DA_p_DB = 0.0;
	double DA_p_DB_r = 0.0;
	double D_R = 0.0;
	double ratio = 0.0;
	double alpha = 0.0;
	double beta = 0.0;
	// double alpha_p_beta = 0.0;
	// double alpha_p_beta_r = 0.0;
	double dx = 0.0;
	double dy = 0.0;
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
	
	p1_index = particle_list[p1].index;
    p2_index = particle_list[p2].index;

    
	D_A = species_list[p1_index].DiffCoeff;
	D_B = species_list[p2_index].DiffCoeff;

	DA_p_DB = D_A + D_B;
	DA_p_DB_r = 1.0/DA_p_DB;

	alpha = D_B*DA_p_DB_r;
	beta = D_A*DA_p_DB_r;
	D_R = D_A*D_B*DA_p_DB_r;

	// alpha = sqrt(D_B/D_A);
	// beta = sqrt(D_A/D_B);
	// D_R = D_A + D_B;


	x_A = particle_list[p1].x;
	y_A = particle_list[p1].y;
	z_A = particle_list[p1].z;

	x_B = particle_list[p2].x;
	y_B = particle_list[p2].y;
	z_B = particle_list[p2].z;

	x_R = alpha*x_A  + beta*x_B;
	y_R = alpha*y_A  + beta*y_B;
	z_R = alpha*z_A  + beta*z_B;




	dx = gsl_ran_gaussian(gen, sqrt(2*D_R*sampled_time));                
    dy = gsl_ran_gaussian(gen, sqrt(2*D_R*sampled_time));
    dz = gsl_ran_gaussian(gen, sqrt(2*D_R*sampled_time));


    //The locaiton of particle C after bimol1, which is located in the position of center of diffusion vector
    x_R =  x_R  + dx;
    y_R =  y_R  + dy;
    z_R =  z_R  + dz;

   
    while (sttoint(species_list[p1_index].asso[index].partner.c_str(), species_list) != p2_index){
		index++;
	}

	particle_list[p1].x = x_R;
	particle_list[p1].y = y_R;
	particle_list[p1].z = z_R;
	particle_list[p1].r = sqrt(particle_list[p1].x*particle_list[p1].x + particle_list[p1].y*particle_list[p1].y + particle_list[p1].z*particle_list[p1].z);
	if(particle_list[p1].x < 0.0){
		particle_list[p1].phi = atan(particle_list[p1].y/particle_list[p1].x) + Pi;	
	}
	else{
		particle_list[p1].phi = atan(particle_list[p1].y/particle_list[p1].x);
	}
    particle_list[p1].theta = acos(particle_list[p1].z/particle_list[p1].r);
	particle_list[p1].index = sttoint(species_list[p1_index].asso[index].product.c_str(), species_list);
	// std::cout << "particle in element " << p1 << " after bimol1 is: " <<  species_list[particle_list[p1].index].name << std::endl;
	particle_list[p1].update_flag ++;
	particle_list[p2].exist_flag = 0;
	// std::cout << "size of particle_list after bimol1: " << particle_list.size() << std::endl;

}

void bimol1_unbound(std::vector<Species>& species_list,  std::vector<Particle>& particle_list, int p1, int p2, double sampled_time, double sampled_r, double sampled_theta, gsl_rng * gen){


	int index = -1;
	int p1_index = -1;
	int p2_index = -1;
	
	double D_A = 0.0;
	double D_B = 0.0;
	double D_C = 0.0;
	double DA_p_DB = 0.0;
	double DA_p_DB_r = 0.0;
	double D_R = 0.0;
	double ratio = 0.0;
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


	double phi = 0.0;
	// double theta = 0.0;
	double x_r = 0.0;
	double y_r = 0.0;
	double z_r = 0.0;
	double r_r = 0.0;
	double r_r_r = 0.0;

	double x_r_init = 0.0;
	double y_r_init = 0.0;
	double z_r_init = 0.0;
	double r_init = 0.0;
	double theta_init = 0.0;

	p1_index = particle_list[p1].index;
    p2_index = particle_list[p2].index;

	D_A = species_list[p1_index].DiffCoeff;
	D_B = species_list[p2_index].DiffCoeff;
	// std::cout << "radii: " << species_list[particle_list[p1].index].radius + species_list[particle_list[p2].index].radius << std::endl;

	phi = 2*Pi*(double)gsl_rng_get(gen)/(double)gsl_rng_max(gen); // azimuthal angle


	
	
	DA_p_DB = D_A + D_B;
	DA_p_DB_r = 1.0/DA_p_DB;
	alpha = D_B*DA_p_DB_r;
	beta = D_A*DA_p_DB_r;
	D_R = D_A*D_B*DA_p_DB_r;

	alpha_p_beta = alpha + beta;
	alpha_p_beta_r = 1.0/alpha_p_beta;

	x_A = particle_list[p1].x;
	y_A = particle_list[p1].y;
	z_A = particle_list[p1].z;

	// std::cout << x_A << SP << y_A << SP << z_A << std::endl;

	x_B = particle_list[p2].x;
	y_B = particle_list[p2].y;
	z_B = particle_list[p2].z;

	// std::cout << x_B << SP << y_B << SP << z_B << std::endl;

	x_R = alpha * x_A  + beta * x_B;
	y_R = alpha * y_A  + beta * y_B;
	z_R = alpha * z_A  + beta * z_B;

	// std::cout << x_R << SP << y_R << SP << z_R << std::endl;
	

	dx = gsl_ran_gaussian(gen, sqrt(2*D_R*sampled_time));                
    dy = gsl_ran_gaussian(gen, sqrt(2*D_R*sampled_time));
    dz = gsl_ran_gaussian(gen, sqrt(2*D_R*sampled_time));


    //The locaiton of center of diffusion
    x_R =  x_R  + dx;
    y_R =  y_R  + dy;
    z_R =  z_R  + dz;

    x_r_init = x_B - x_A;
    y_r_init = y_B - y_A;
    z_r_init = z_B - z_A;
    r_init = sqrt(x_r_init * x_r_init + y_r_init * y_r_init + z_r_init * z_r_init);


    theta_init = acos(z_r_init / r_init);
    double phi_init = atan(y_r_init / x_r_init);
    
    if (x_r_init < 0.0){
    	phi_init = Pi + phi_init;
    }


    double x_s = sampled_r * sin(sampled_theta) * cos(phi); 
    double y_s = sampled_r * sin(sampled_theta) * sin(phi); 
    double z_s = sampled_r * cos(sampled_theta); 



    //Y rotation;
    double x_sss = x_s * cos(theta_init) + z_s * sin(theta_init);
    double y_sss = y_s;
    double z_sss = -x_s * sin(theta_init) + z_s * cos(theta_init);


    //Z rotation
    x_r = x_sss * cos(phi_init) - y_sss * sin(phi_init);
    y_r = x_sss * sin(phi_init) + y_sss * cos(phi_init);
    z_r = z_sss;



// //X rotation

//     x_r = x_r;
//     y_r = y_r * cos(Beta) - z_r * sin(Beta);
//     z_r = y_r * sin(Beta) + z_r * cos(Beta);

	particle_list[p1].x = (x_R - beta * x_r) * alpha_p_beta_r;
	particle_list[p1].y = (y_R - beta * y_r) * alpha_p_beta_r;
	particle_list[p1].z = (z_R - beta * z_r) * alpha_p_beta_r;
	particle_list[p1].r = sqrt(particle_list[p1].x * particle_list[p1].x + particle_list[p1].y * particle_list[p1].y + particle_list[p1].z * particle_list[p1].z);
	if(particle_list[p1].x < 0.0){
		particle_list[p1].phi = atan(particle_list[p1].y/particle_list[p1].x) + Pi;	
	}
	else{
		particle_list[p1].phi = atan(particle_list[p1].y/particle_list[p1].x);
	}
	particle_list[p1].theta = acos(particle_list[p1].z/particle_list[p1].r);

	// std::cout << particle_list[p1].x << SP << particle_list[p1].y << SP << particle_list[p1].z << std::endl;

	particle_list[p2].x = (x_R + alpha * x_r) * alpha_p_beta_r;
	particle_list[p2].y = (y_R + alpha * y_r) * alpha_p_beta_r;
	particle_list[p2].z = (z_R + alpha * z_r) * alpha_p_beta_r;

	particle_list[p2].r = sqrt(particle_list[p2].x*particle_list[p2].x + particle_list[p2].y*particle_list[p2].y + particle_list[p2].z*particle_list[p2].z);
	if(particle_list[p2].x < 0.0){
		particle_list[p2].phi = atan(particle_list[p2].y/particle_list[p2].x) + Pi;	
	}
	else{
		particle_list[p2].phi = atan(particle_list[p2].y/particle_list[p2].x);
	}
	
	particle_list[p2].theta = acos(particle_list[p2].z/particle_list[p2].r);

	// std::cout << particle_list[p2].x << SP << particle_list[p2].y << SP << particle_list[p2].z << std::endl;
	particle_list[p1].update_flag ++;
	particle_list[p2].update_flag ++;


}

