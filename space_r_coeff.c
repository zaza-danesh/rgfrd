
//
//  space_r_coeff.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/05/05.
//
//

#include "space_r_coeff.h"


std::vector<space_r_coeff_struct> space_r_coeff(std::vector<Species>& species_list, std::vector<Roots>& species_roots, double b_value ){


	
	double ka = 0.0; double D = 0.0; double sigma = 0.0; double b = 0.0; double kd = 0.0; double gamma_b = 0.0; double gamma_u = 0.0;
	int J = 0;
	int j = 0;
	int i = 0;
	int ii = 0;
    int index = -1;
	int max_num_root = MAX_NUM_ROOT;
	struct phys_parameters_bimol2 physparams;
	struct space_r_coeff_struct space_r_coeff;
	std::vector<space_r_coeff_struct> space_r_coeff_arr;

	
    double q = 0.0;
    double q_2pi = 0.0;
    double q_2pi_r = 0.0;
    double q_2pi_sqrd = 0.0;
    double q_2pi_kappa = 0.0;
    double q_2pi_kappa_r = 0.0;
    double q_2pi_kappa_sqrd = 0.0;
    double q_2pi_over_b_m_sigma = 0.0;
    double q_2pi_over_b_m_sigma_r = 0.0;
     
    double l = 0.0;
    double l_2pi = 0.0;
    double l_2pi_r = 0.0;
    double l_2pi_sqrd = 0.0;
    double l_2pi_kappa = 0.0;
    double l_2pi_kappa_r = 0.0;
    double l_2pi_kappa_sqrd = 0.0;
    double l_2pi_over_b_m_sigma = 0.0;
    double l_2pi_over_b_m_sigma_r = 0.0; 

    double term1_q = 0.0;
    double term2_q = 0.0;
    double term3_q = 0.0;
    double term4_q = 0.0;

    double term1_l = 0.0;
    double term2_l = 0.0;
    double term3_l = 0.0;
    double term4_l = 0.0;

    double psi_y = 0.0;  
    double hd_y = 0.0;  
    double psi_x = 0.0;
    double hd_x = 0.0;  
    double cosine_term = 0.0;
    double cosine_term1 = 0.0;
    double cosine_term2 = 0.0;
    double sine_term = 0.0;
    double sine_term1 = 0.0;
    double sine_term2 = 0.0;
    double coshype_term = 0.0;
    double coshype_term1 = 0.0;
    double coshype_term2 = 0.0;
    double sinhype_term = 0.0;
    double sinhype_term1 = 0.0;
    double sinhype_term2 = 0.0;
    double denominator_y = 0.0;
    double denominator_x = 0.0;
    double exp_term_y = 0.0;
    double exp_term_x = 0.0;
    double partial_sum_y = 0.0;
    double sum_y = 0.0;
    double f_x = 0.0;

    double result = 0.0;


    double sqrt_D = 0.0;
    double sqrt_D_r = 0.0;
    double b_m_sigma = 0.0;
    double b_m_sigma_r = 0.0;
    double kappa = 0.0;
    double kappa1 = 0.0;
    double kappa1_r = 0.0;
    double kappa_r = 0.0;
    double kappa_sqrd = 0.0;
    double kD = 0.0;
    double kD_r = 0.0;
    double N_disso = 0.0;

	for (i = 0 ; i < species_roots.size() ; i++){
            index = species_roots[i].index;
			space_r_coeff.index = index;			
			ka = species_list[index].disso[J].disso_rate.ka;
			D = species_list[sttoint(species_list[index].disso[J].product1.c_str(), species_list)].DiffCoeff +species_list[sttoint(species_list[index].disso[J].product2.c_str(), species_list)].DiffCoeff;
			sigma = species_list[sttoint(species_list[index].disso[J].product1.c_str(), species_list)].radius + species_list[sttoint(species_list[index].disso[J].product2.c_str(), species_list)].radius;
			// b = std::min( ALPHA*sigma, b_value + sigma);
            b = ALPHA*sigma;
            space_r_coeff.b = b;
			kd = species_list[index].disso[J].disso_rate.kd;
			gamma_b = species_list[index].disso[J].disso_rate.gamma_b;
			gamma_u = species_list[index].disso[J].disso_rate.gamma_u;
			sqrt_D = sqrt(D);
			sqrt_D_r = 1.0/sqrt_D;
			b_m_sigma = b - sigma;
			b_m_sigma_r = 1./b_m_sigma;
			kappa = sqrt_D*b_m_sigma_r;
			kappa1 = sqrt_D/sigma;    
			kappa1_r = 1./kappa1;
			kappa_r = 1./kappa;
			kappa_sqrd = kappa*kappa;
			kD = 4.0*Pi*sigma*D;
			kD_r = 1.0/kD;
			N_disso = 1.0 + ka*kD_r;   
			space_r_coeff.pre_factor = -kd/D;


			if (species_roots[i].real_root != 0){
    		    max_num_root--;
                space_r_coeff.real_root_flag = 1;
    		    l = species_roots[i].real_root;
    		    l_2pi = 2.0*Pi*l;
    		    l_2pi_kappa = kappa * l_2pi;
    		    l_2pi_kappa_r = 1./l_2pi_kappa;
    		    l_2pi_kappa_sqrd = l_2pi_kappa*l_2pi_kappa;
    		    l_2pi_r = 1.0/l_2pi;
    		    l_2pi_over_b_m_sigma = l_2pi*b_m_sigma_r;
    		    l_2pi_over_b_m_sigma_r = 1./l_2pi_over_b_m_sigma;
    		    l_2pi_sqrd = l_2pi*l_2pi;

                term1_l = l_2pi_kappa_sqrd - gamma_u + gamma_b;
                term1_l = term1_l * l_2pi_kappa_r;
                term1_l = term1_l * (kappa1_r + N_disso*kappa_r);
                term1_l = term1_l + 2*l_2pi_kappa * kappa1_r;
                term1_l = term1_l + kd * (kappa1_r + kappa_r)*l_2pi_kappa_r;
                term2_l = l_2pi_kappa_sqrd - gamma_u + gamma_b;
                term2_l = N_disso * term2_l;
                term2_l = term2_l + kd;
                term3_l = l_2pi_kappa_sqrd -gamma_u + gamma_b + kd;
                term3_l = term3_l * kappa1_r * l_2pi_kappa;
                term3_l = 1./term3_l;
                term4_l = l_2pi_kappa_sqrd -gamma_u + gamma_b + kd;
                term4_l = term4_l * kappa_r * kappa1_r;
                term4_l = term4_l + 2*N_disso;
                psi_x = term1_l * term2_l * term3_l  - term4_l;
                psi_x = -0.5*psi_x;
                coshype_term = cosh(l_2pi);
                sinhype_term = sinh(l_2pi);
                hd_x = psi_x*sinhype_term;
                space_r_coeff.denominator_x = 1.0/hd_x;
                space_r_coeff.coshype_p_sinhype = sigma * l_2pi_over_b_m_sigma_r * coshype_term + l_2pi_over_b_m_sigma_r * l_2pi_over_b_m_sigma_r * sinhype_term;
                space_r_coeff.coshype_sinhype_arg = l_2pi_over_b_m_sigma;    
                space_r_coeff.coshype_sinhype_coeff = l_2pi_over_b_m_sigma_r;    
                space_r_coeff.exp_arg_x = gamma_u - l_2pi_kappa_sqrd;

				ii++;
    		}
    		else{
    		    ii = 0;
    		}
    		for (j = 0 ; j < max_num_root ; j++){
				q = species_roots[i].char_roots[j];
        		q_2pi = 2.0*Pi*q;
        		q_2pi_kappa = kappa*q_2pi;
        		q_2pi_kappa_r = 1./q_2pi_kappa;
        		q_2pi_kappa_sqrd = q_2pi_kappa*q_2pi_kappa;
        		q_2pi_r = 1.0/q_2pi;
        		q_2pi_over_b_m_sigma = q_2pi*b_m_sigma_r;
        		q_2pi_over_b_m_sigma_r = 1./q_2pi_over_b_m_sigma;
        		q_2pi_sqrd = q_2pi*q_2pi;        
        		term1_q = -q_2pi_kappa_sqrd - gamma_u + gamma_b;
                term1_q = term1_q * q_2pi_kappa_r;
                term1_q = term1_q * (kappa1_r + N_disso*kappa_r);
                term1_q = term1_q - 2*q_2pi_kappa * kappa1_r;
                term1_q = term1_q + kd * (kappa1_r + kappa_r)*q_2pi_kappa_r;
                term2_q = -q_2pi_kappa_sqrd - gamma_u + gamma_b;
                term2_q = N_disso * term2_q;
                term2_q = term2_q + kd;
                term3_q = -q_2pi_kappa_sqrd -gamma_u + gamma_b + kd;
                term3_q = term3_q * kappa1_r*q_2pi_kappa;
                term3_q = 1./term3_q;
                term4_q= -q_2pi_kappa_sqrd -gamma_u + gamma_b + kd;
                term4_q = term4_q * kappa_r*kappa1_r;
                term4_q = term4_q + 2*N_disso;
                psi_y = term1_q * term2_q * term3_q  + term4_q;
                psi_y = 0.5*psi_y;
                cosine_term = cos(q_2pi);
                sine_term = sin(q_2pi);
                hd_y = psi_y*sine_term;        
                space_r_coeff.denominator_y[ii] = 1.0/hd_y;
                space_r_coeff.cosine_p_sine[ii] = sigma * q_2pi_over_b_m_sigma_r * cosine_term + q_2pi_over_b_m_sigma_r * q_2pi_over_b_m_sigma_r * sine_term ;
                space_r_coeff.cosine_sine_arg[ii] = q_2pi_over_b_m_sigma ; 
                space_r_coeff.cosine_sine_coeff[ii] = q_2pi_over_b_m_sigma_r ; 
                space_r_coeff.exp_arg_y[ii] = gamma_u + q_2pi_kappa_sqrd;
        		ii++;
        	}
        space_r_coeff_arr.push_back(space_r_coeff);

	}

	return space_r_coeff_arr;
}

