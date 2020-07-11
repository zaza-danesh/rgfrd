
//
//  Qt_coeff.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/05/05.
//
//

#include "Qt_coeff.h"


std::vector<Qt_coeff_struct> Qt_coeff(std::vector<Species>& species_list, std::vector<Roots>& species_roots,  double b_value ){


	
	double ka = 0.0; double D = 0.0; double sigma = 0.0; double b = 0.0; double kd = 0.0; double gamma_b = 0.0; double gamma_u = 0.0;
	int J = 0;
    int j = 0;
	int i = 0;
	int ii = 0;
    int index = -1;
	int max_num_root = MAX_NUM_ROOT;
	struct phys_parameters_bimol2 physparams;
	struct Qt_coeff_struct Qt_coeff;
	std::vector<Qt_coeff_struct> Qt_coeff_arr;
	
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
    double Term1_q = 0.0;
    double term2_q = 0.0;
    double Term2_q = 0.0;
    double term3_q = 0.0;
    double term4_q = 0.0;

    double term1_l = 0.0;
    double Term1_l = 0.0;
    double term2_l = 0.0;
    double Term2_l = 0.0;
    double term3_l = 0.0;
    double term4_l = 0.0;

    double psi_y = 0.0;  
    double Psi_y = 0.0;
    double hd_y = 0.0;  
    double psi_x = 0.0;
    double Psi_x = 0.0;
    double hd_x = 0.0;  
    double cosine_term = 0.0;
    double sine_term = 0.0;
    double coshype_term = 0.0;
    double sinhype_term = 0.0;
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
			Qt_coeff.index = index;			
			ka = species_list[index].disso[J].disso_rate.ka;
			D = species_list[sttoint(species_list[index].disso[J].product1.c_str(), species_list)].DiffCoeff +species_list[sttoint(species_list[index].disso[J].product2.c_str(), species_list)].DiffCoeff;
			sigma = species_list[sttoint(species_list[index].disso[J].product1.c_str(), species_list)].radius + species_list[sttoint(species_list[index].disso[J].product2.c_str(), species_list)].radius;
			// b = std::min( ALPHA*sigma, b_value + sigma);
            b = ALPHA*sigma;
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
			Qt_coeff.pre_factor1 = kd;
            Qt_coeff.pre_factor2 =  gamma_b;
			// physparams = make_parameters_bimol2(ka, D, sigma, b , kd, gamma_b, gamma_u );
			// find_roots(roots.real_root, roots.char_roots,  physparams);

			if (species_roots[i].real_root != 0){
    		    max_num_root--;
    		    l = species_roots[i].real_root;
                l_2pi = 2.0*Pi*l;
                l_2pi_kappa = kappa * l_2pi;
                l_2pi_kappa_r = 1./l_2pi_kappa;
                l_2pi_kappa_sqrd = l_2pi_kappa*l_2pi_kappa;
                l_2pi_r = 1.0/l_2pi;
                l_2pi_over_b_m_sigma = l_2pi*b_m_sigma_r;
                l_2pi_over_b_m_sigma_r = 1./l_2pi_over_b_m_sigma;
                l_2pi_sqrd = l_2pi*l_2pi;
                coshype_term = cosh(l_2pi);
                sinhype_term = sinh(l_2pi);
                term1_l = l_2pi_kappa_sqrd - gamma_u + gamma_b;
                term1_l = term1_l * l_2pi_kappa_r;
                term1_l = term1_l * (kappa1_r + N_disso*kappa_r);
                term1_l = term1_l + 2.0*l_2pi_kappa * kappa1_r;
                term1_l = term1_l + kd * (kappa1_r + kappa_r)*l_2pi_kappa_r;
                Term1_l = term1_l*coshype_term;
                term4_l = l_2pi_kappa_sqrd -gamma_u + gamma_b + kd;
                term4_l = term4_l * kappa_r * kappa1_r;
                term4_l = term4_l + 2.0*N_disso;
                Term2_l = term4_l*sinhype_term;
                Psi_x = Term1_l + Term2_l;
                Psi_x = 0.5*Psi_x;
                hd_x = Psi_x;        
                denominator_x = gamma_u - gamma_b - l_2pi_kappa_sqrd;
                denominator_x = denominator_x * hd_x;
                denominator_x = 1./denominator_x;
                coshype_term = sigma * l_2pi_over_b_m_sigma * coshype_term;
    		    Qt_coeff.lambda1[ii] = (gamma_u - gamma_b - l_2pi_kappa_sqrd);
    		    Qt_coeff.lambda2[ii] = denominator_x * (coshype_term + sinhype_term);				
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
                cosine_term = cos(q_2pi);
                sine_term = sin(q_2pi);
                term1_q = -q_2pi_kappa_sqrd - gamma_u + gamma_b;
                term1_q = term1_q * q_2pi_kappa_r;
                term1_q = term1_q * (kappa1_r + N_disso*kappa_r);
                term1_q = term1_q - 2.0*q_2pi_kappa * kappa1_r;
                term1_q = term1_q + kd * (kappa1_r + kappa_r)*q_2pi_kappa_r;
                Term1_q = -term1_q*cosine_term;
                term4_q= -q_2pi_kappa_sqrd -gamma_u + gamma_b + kd;
                term4_q = term4_q * kappa_r*kappa1_r;
                term4_q = term4_q + 2.0*N_disso;
                Term2_q = term4_q*sine_term;
                Psi_y = Term1_q + Term2_q;
                Psi_y = 0.5*Psi_y;
                hd_y = Psi_y;
                denominator_y = gamma_u - gamma_b + q_2pi_kappa_sqrd;
                denominator_y = denominator_y*hd_y;
                denominator_y = 1./denominator_y;    
                cosine_term = sigma * q_2pi_over_b_m_sigma * cosine_term;
        		Qt_coeff.lambda1[ii] = (gamma_u - gamma_b + q_2pi_kappa_sqrd);
        		Qt_coeff.lambda2[ii] = denominator_y * (cosine_term + sine_term);
        		ii++;
        	}
        Qt_coeff_arr.push_back(Qt_coeff);
	}

	return Qt_coeff_arr;
}

