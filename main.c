//
//  main_rgfrd.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/21.
//
//
 

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
int KKKK = 0;

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
	double final_time = 1;
	double sampled_time = 0.0;
	double t_previous = 0.0;
	double b_value = 0.0;
	double max_diffusion_time_step = 0.0;
	double seed = atof(argv[1]);

	double TT = atof(argv[2]);
	double H = 2.0;
	double R_curv = 5;
	char negative_distances = 'T';
	char three_body_problem = 'F';
	double ref_flag = atof(argv[11]);
	std::string event = "--";

	gsl_rng * gen;  //generator
    const gsl_rng_type * T;
    // gsl_rng_env_setup(); 
    // final_time = atof(argv[1]);
    T = gsl_rng_default; //default is mt19937
    gen = gsl_rng_alloc (T);
    gsl_rng_set(gen, seed);


    double S_init = atoi(argv[3])/2;
	double P_init = atoi(argv[3])/2;
	double E1_init = atof(argv[4]);
	double E2_init = atof(argv[4]);
	double SE1_init = 0;
	double PE2_init = 0;
	double SP_tot = S_init + P_init;


	double sigma = atof(argv[5]);
	double Rp = 2*sigma; //2sigma
    double diffusion_coefficient =  atof(argv[6]);
    double Dp = 2*diffusion_coefficient; //2D
    double kD = 4*Pi*Rp*Dp;

    double sphere_radius = atof(argv[10]);
    double Vol = 4.0*Pi*(sphere_radius*sphere_radius*sphere_radius)/3.0;


    double ka1 = atof(argv[7]);
    std::cout <<"Nd: " <<  1 + ka1/kD << std::endl;
    // double kd1 = atof(argv[8]);
    // double k1 = atof(argv[9]);
    // double f = atof(argv[9]);
    // double ka2 = atof(argv[10]);
    double ka2 = ka1;
    double kd2 = atof(argv[8]);
    double k2 = atof(argv[9]);
    double f = atof(argv[12]);
    double k1 = (2.0*f/5.0)*k2;
    double kd1 = (kd2+k2) - k1;
    
    double kon1 = 1.0/(1.0/ka1 + 1.0/kD );
    double kon2 = 1.0/(1.0/ka2 + 1.0/kD );
    double K1 = Vol*(kd1 + k1 )/(kon1*(S_init + P_init));
    double K2 = K1;

    // std::cout << "K: " << K1 << std::endl;
    // std::cout << "K2: " << K2 << std::endl;

if (f==0){
	double delta_x = 0.01;
	double x = 0.0;
	double temp1=0.0;
	double temp2=0.0;
	std::ofstream p_theory; p_theory.open(DATA_OUTPUT_PATH "data_theory.dat");
    while(x < 2 ){
    		temp1 = x-1;
    		temp2 = K2*(K1/K2 + x);
    		p_theory << (temp1 - temp2 + sqrt((temp1 - temp2)*(temp1 - temp2) + 4*K2*temp1*x))/(2*temp1) << std::endl;

    	x = x + delta_x;
    }
    p_theory.close();

}

    
	// keff = 1/(1/kf_K + 1/(4*PI*Dtot*SIGMA));    //Morelli
 //    K1 = (kd_K+kprod_K) / (keff*NStot0/V);		//Morelli
    

   
	double Number_of_particles = S_init + P_init + E1_init + E2_init + SE1_init + PE2_init;
	np = Number_of_particles;

	// std::ofstream system_file;
	// system_file.open(DATA_OUTPUT_PATH "system.txt");
	// system_file << "//System" << std::endl;
	// system_file << "Number_of_particles: " << Number_of_particles  << std::endl;
	// system_file << "Sphere_radius: " << 1.0 << std::endl;
	// system_file.close();
	// std::cout << __LINE__ << std::endl;

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

    // get_system(sphere_radius, np);


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

	get_reactions(get_number_of_species(), species_list, ka1, kd1, k1, ka2, kd2, k2);

	double sum = 0.0;
	// //DEBUGK
	for (i = 1 ; i < get_number_of_species() ; i++){
		sum = sum + species_list[i].quantity * species_list[i].radius * species_list[i].radius * species_list[i].radius;
	}

	particle_list = initial_distribution(species_list, dij_matrix, np, sphere_radius, gen, k);




	// std::ofstream  number_of_particles_file, t_file;
 //    std::ofstream nt_file; nt_file.open("nt.dat");
 //    std::ofstream n2t_file; n2t_file.open("n2t.dat");
 //    std::ofstream n_file;  n_file.open("n.dat");

std::ofstream pt_file; pt_file.open("pt.dat");


    int j = 0;
    int NP = 0;

    	species_roots = Find_Roots(species_list, b_value);
    	It_coeff_arr =	It_coeff(species_list, species_roots, b_value);
    	Qt_coeff_arr =	Qt_coeff(species_list, species_roots, b_value);
    	Jt_coeff_arr =	Jt_coeff(species_list, species_roots, b_value);
    	space_r_coeff_arr =	space_r_coeff(species_list, species_roots, b_value);
	

double S = 0.0;
double P = 0.0;
double E1 = 0.0;
double E2 = 0.0;
double SE1 = 0.0;
double PE2 = 0.0;


double St = 0.0;
double Pt = 0.0;
double E1t = 0.0;
double E2t = 0.0;
double SE1t = 0.0;
double PE2t = 0.0;

double S2t = 0.0;
double P2t = 0.0;
double E12t = 0.0;
double E22t = 0.0;
double SE12t = 0.0;
double PE22t = 0.0;


// double Ct = 0.0;
// double A2t = 0.0;
// double B2t = 0.0;
// double C2t = 0.0;
double t1 = 0.0;
int N_threebodyproblem =0;
double tt = 0;
double time_to_save = 5e6;

	// std::cout << __LINE__ << std::endl;


    // while( 0 < 1){
    while( ref_flag > 0){
		k++;
    	if ( k != 0){
    		if (np != particle_list.size()){
    			np = particle_list.size();
    			dij_matrix.resize(np);
    			for (i = 0 ; i < dij_matrix.size() ; i++){
    				dij_matrix[i].resize(np);	
    			}    			
    		}		
    		negative_distances = calculate_distances(species_list, particle_list, dij_matrix, np, k);
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
    	particle_particle_diffusion_time_array = particle_particle_diffusion(species_list, particle_list, dij_matrix, H, k );
    	particle_boundary_diffusion_time_array = particle_boundary_diffusion(species_list, particle_list, dij_matrix, H, sphere_radius, R_curv, k );
    	pair_list = make_pairs(three_body_problem, number_of_pairs, max_diffusion_time_step, b_value, species_list, particle_list,  particle_particle_diffusion_time_array, particle_boundary_diffusion_time_array[0].t_w, sphere_radius, H, k );
    	// std::cout << number_of_pairs << std::endl;
    	if (three_body_problem == 'T'){
    		// N_threebodyproblem++;
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
						if(t > time_to_save){
							t1 = t1 + sampled_time;
							tt = tt + sampled_time;
							S=0; P=0; E1=0; E2=0; SE1=0; PE2=0;
							for (i = 0 ; i < particle_list.size() ; i++){
								if (particle_list[i].index == 1){
									S++;
								}
								else if(particle_list[i].index == 2){
									P++;
								}
								else if(particle_list[i].index == 3){
									E1++;
								}
								else if(particle_list[i].index == 4){
									E2++;
								}
								else if(particle_list[i].index == 5){
									SE1++;
								}
								else if(particle_list[i].index == 6){
									PE2++;
								}
							}
						}
						// std::cout << A << SP << B << SP << C << std::endl;
						St = St + S * sampled_time;
						Pt = Pt + P * sampled_time;
						E1t = E1t + E1 * sampled_time;
						E2t = E2t + E2 * sampled_time;
						SE1t = SE1t + SE1 * sampled_time;
						PE2t = PE2t + PE2 * sampled_time;
						S2t = S2t + S * S * sampled_time;
						P2t = P2t + P * P * sampled_time;
						E12t = E12t + E1 * E1 * sampled_time;
						E22t = E22t + E2 * E2 * sampled_time;
						SE12t = SE12t + SE1 * SE1 * sampled_time;
						PE22t = PE22t + PE2 * PE2 * sampled_time;
						update_system_bimolecular_reaction_type_one(species_list, particle_list, pair_list, reaction_event.pair_index,  It_coeff_arr, Qt_coeff_arr, space_r_coeff_arr,  b_value, sampled_time, gen);
						reflection(species_list, particle_list, sphere_radius, ref_flag, final_time, k);
						if(t1 > TT){
							t1 = 0.0;
							pt_file << (Pt/tt + PE2t/tt)/SP_tot << std::endl;
							// n_file << S << SP << P << SP << E1 << SP << E2 << SP << SE1 << SP << PE2 << std::endl;
							// nt_file << St/tt << SP << Pt/tt <<  SP << E1t/tt << SP << E2t/tt << SP << SE1t/tt << SP << PE2t/tt << std::endl;
							// n2t_file << S2t/tt - (St/tt)*(St/tt) << SP << P2t/tt -(Pt/tt)*(Pt/tt) <<  SP << E12t/tt - (E1t/tt)*(E1t/tt) << SP << E22t/tt - (E2t/tt)*(E2t/tt) << SP << SE12t/tt - (SE1t/tt)*(SE1t/tt) << SP << PE22t/tt - (PE2t/tt)*(PE2t/tt) << std::endl;
						}
					}//else if bimol1
					else if ( event == "bimol2" ){
						t = t + sampled_time;						
						if(t > time_to_save){
							t1 = t1 + sampled_time;
							tt = tt + sampled_time;
							S=0; P=0; E1=0; E2=0; SE1=0; PE2=0;
							for (i = 0 ; i < particle_list.size() ; i++){
								if (particle_list[i].index == 1){
									S++;
								}
								else if(particle_list[i].index == 2){
									P++;
								}
								else if(particle_list[i].index == 3){
									E1++;
								}
								else if(particle_list[i].index == 4){
									E2++;
								}
								else if(particle_list[i].index == 5){
									SE1++;
								}
								else if(particle_list[i].index == 6){
									PE2++;
								}
							}
						}
						St = St + S * sampled_time;
						Pt = Pt + P * sampled_time;
						E1t = E1t + E1 * sampled_time;
						E2t = E2t + E2 * sampled_time;
						SE1t = SE1t + SE1 * sampled_time;
						PE2t = PE2t + PE2 * sampled_time;
						S2t = S2t + S * S * sampled_time;
						P2t = P2t + P * P * sampled_time;
						E12t = E12t + E1 * E1 * sampled_time;
						E22t = E22t + E2 * E2 * sampled_time;
						SE12t = SE12t + SE1 * SE1 * sampled_time;
						PE22t = PE22t + PE2 * PE2 * sampled_time;
						update_system_bimolecular_reaction_type_two(k, species_list, particle_list, pair_list, reaction_event.pi_bimol2,  It_coeff_arr, Qt_coeff_arr, Jt_coeff_arr, space_r_coeff_arr,  b_value, sampled_time, gen);
						reflection(species_list, particle_list, sphere_radius, ref_flag,  final_time, k);
						if(t1 > TT){
							t1 = 0.0;
							pt_file << (Pt/tt + PE2t/tt)/SP_tot << std::endl;
							// n_file << S << SP << P << SP << E1 << SP << E2 << SP << SE1 << SP << PE2 << std::endl;
							// nt_file << St/tt << SP << Pt/tt <<  SP << E1t/tt << SP << E2t/tt << SP << SE1t/tt << SP << PE2t/tt << std::endl;
							// n2t_file << S2t/tt - (St/tt)*(St/tt) << SP << P2t/tt -(Pt/tt)*(Pt/tt) <<  SP << E12t/tt - (E1t/tt)*(E1t/tt) << SP << E22t/tt - (E2t/tt)*(E2t/tt) << SP << SE12t/tt - (SE1t/tt)*(SE1t/tt) << SP << PE22t/tt - (PE2t/tt)*(PE2t/tt) << std::endl;
						}
					}//elsif bimol1
					else if ( event == "monomol" ){
						t = t + sampled_time;
						if(t > time_to_save){
							t1 = t1 + sampled_time;
							tt = tt + sampled_time;
							S=0; P=0; E1=0; E2=0; SE1=0; PE2=0;
							for (i = 0 ; i < particle_list.size() ; i++){
								if (particle_list[i].index == 1){
									S++;
								}
								else if(particle_list[i].index == 2){
									P++;
								}
								else if(particle_list[i].index == 3){
									E1++;
								}
								else if(particle_list[i].index == 4){
									E2++;
								}
								else if(particle_list[i].index == 5){
									SE1++;
								}
								else if(particle_list[i].index == 6){
									PE2++;
								}
							}
						}
						// std::cout << A << SP << B << SP << C << std::endl;
						St = St + S * sampled_time;
						Pt = Pt + P * sampled_time;
						E1t = E1t + E1 * sampled_time;
						E2t = E2t + E2 * sampled_time;
						SE1t = SE1t + SE1 * sampled_time;
						PE2t = PE2t + PE2 * sampled_time;
						S2t = S2t + S * S * sampled_time;
						P2t = P2t + P * P * sampled_time;
						E12t = E12t + E1 * E1 * sampled_time;
						E22t = E22t + E2 * E2 * sampled_time;
						SE12t = SE12t + SE1 * SE1 * sampled_time;
						PE22t = PE22t + PE2 * PE2 * sampled_time;
						update_system_monomolecular_reaction(species_list, particle_list, pair_list, reaction_event.pi_monomol,  It_coeff_arr, Qt_coeff_arr, space_r_coeff_arr,  b_value, sampled_time, gen);
						reflection(species_list, particle_list, sphere_radius, ref_flag,  final_time, k);				
						if(t1 > TT){
							t1 = 0.0;
							pt_file << (Pt/tt + PE2t/tt)/SP_tot << std::endl;
							// n_file << S << SP << P << SP << E1 << SP << E2 << SP << SE1 << SP << PE2 << std::endl;
							// nt_file << St/tt << SP << Pt/tt <<  SP << E1t/tt << SP << E2t/tt << SP << SE1t/tt << SP << PE2t/tt << std::endl;
							// n2t_file << S2t/tt - (St/tt)*(St/tt) << SP << P2t/tt -(Pt/tt)*(Pt/tt) <<  SP << E12t/tt - (E1t/tt)*(E1t/tt) << SP << E22t/tt - (E2t/tt)*(E2t/tt) << SP << SE12t/tt - (SE1t/tt)*(SE1t/tt) << SP << PE22t/tt - (PE2t/tt)*(PE2t/tt) << std::endl;
						}
					} //elsif monomol
					else {
						std::cout << __FILE__ << SP << __LINE__ <<  ": Error: no event chosen." << std::endl;	
					}
				}
				else {
					t = t + sampled_time;
					if(t > time_to_save){
						t1 = t1 + sampled_time;
						tt = tt + sampled_time;
						S=0; P=0; E1=0; E2=0; SE1=0; PE2=0;
						for (i = 0 ; i < particle_list.size() ; i++){
							if (particle_list[i].index == 1){
								S++;
							}
							else if(particle_list[i].index == 2){
								P++;
							}
							else if(particle_list[i].index == 3){
								E1++;
							}
							else if(particle_list[i].index == 4){
								E2++;
							}
							else if(particle_list[i].index == 5){
								SE1++;
							}
							else if(particle_list[i].index == 6){
								PE2++;
							}
						}
						}
					// std::cout << A << SP << B << SP << C << std::endl;
					St = St + S * sampled_time;
					Pt = Pt + P * sampled_time;
					E1t = E1t + E1 * sampled_time;
					E2t = E2t + E2 * sampled_time;
					SE1t = SE1t + SE1 * sampled_time;
					PE2t = PE2t + PE2 * sampled_time;
					S2t = S2t + S * S * sampled_time;
					P2t = P2t + P * P * sampled_time;
					E12t = E12t + E1 * E1 * sampled_time;
					E22t = E22t + E2 * E2 * sampled_time;
					SE12t = SE12t + SE1 * SE1 * sampled_time;
					PE22t = PE22t + PE2 * PE2 * sampled_time;
					update_system_no_reaction(k, species_list, particle_list, pair_list,  It_coeff_arr, Qt_coeff_arr, space_r_coeff_arr,  b_value, sampled_time, gen);
					reflection(species_list, particle_list, sphere_radius, ref_flag,  final_time, k);
					if(t1 > TT){
						t1 = 0.0;
						pt_file << (Pt/tt + PE2t/tt)/SP_tot << std::endl;
						// n_file << S << SP << P << SP << E1 << SP << E2 << SP << SE1 << SP << PE2 << std::endl;
						// nt_file << St/tt << SP << Pt/tt <<  SP << E1t/tt << SP << E2t/tt << SP << SE1t/tt << SP << PE2t/tt << std::endl;
						// n2t_file << S2t/tt - (St/tt)*(St/tt) << SP << P2t/tt -(Pt/tt)*(Pt/tt) <<  SP << E12t/tt - (E1t/tt)*(E1t/tt) << SP << E22t/tt - (E2t/tt)*(E2t/tt) << SP << SE12t/tt - (SE1t/tt)*(SE1t/tt) << SP << PE22t/tt - (PE2t/tt)*(PE2t/tt) << std::endl;
					}
				}//else, no reaction happens
			}//if check_max_time
    		else{
    			t = t + final_time;
    		}
    	}
    }//while
    
    pt_file.close();
    // nt_file.close();
    // n2t_file.close();
    // n_file.close();

}//main
