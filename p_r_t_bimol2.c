//
//  p_r_t_bimol2.h
//  
//
//  Created by Zahedeh Bashardanesh on 05/27/14.
//
//

 

#include "p_r_t_bimol2.h"



double p_r_t_bimol2(double r, double t, double ka, double D, double sigma, double b, double kd, double gamma_b, double gamma_u, double real_root, double* char_roots ){

	int i = 0;
    int max_num_root = MAX_NUM_ROOT;

    struct phys_parameters_bimol2 physparams;
        physparams.ka = ka;
        physparams.D = D;
        physparams.sigma = sigma;
        physparams.b = b;
        physparams.kd = kd;
        physparams.gamma_b = gamma_b;
        physparams.gamma_u = gamma_u;

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
    double sine_term1 = 0.0;
    double sine_term2 = 0.0;
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

    const double sqrt_D = sqrt(D);
    const double sqrt_D_r = 1.0/sqrt_D;
    const double b_m_sigma = b - sigma;
    const double b_m_sigma_r = 1./b_m_sigma;
    const double kappa1_r = sigma*sqrt_D_r; 
    const double kappa2_r = b*sqrt_D_r;
    const double kappa_r = kappa2_r - kappa1_r;
    const double kappa1 = 1./kappa1_r;
    const double kappa2 = 1./kappa2_r;
    const double kappa = 1./kappa_r;
    const double kappa_sqrd = kappa*kappa;
    const double kD = 4*Pi*sigma*D;
    const double kD_r = 1.0/kD;
    const double N_disso = 1 + ka*kD_r;
    const double factor = -kd*sigma*kD_r;
    const double r_r = 1./r;
    
     if (real_root !=0){
        max_num_root--;
    }


    for (i=0; i<max_num_root ; i++){

        q = char_roots[i];
        q_2pi = 2*Pi*q;
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
        
  
        sine_term1 = sin(q_2pi);

        hd_y = psi_y*sine_term1;
        denominator_y = 1.0/hd_y;



        exp_term_y = exp(-(gamma_u + q_2pi_kappa_sqrd)*t);

        
        sine_term2 =  sin(q_2pi_over_b_m_sigma * (r-b))*r_r;

        partial_sum_y = exp_term_y * sine_term2 * denominator_y;
      

        sum_y = sum_y + partial_sum_y;
        
              
    }

    if (real_root != 0){

        l = real_root;
        l_2pi = 2*Pi*l;
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

        
        sinhype_term1 = sinh(l_2pi);

        hd_x = psi_x*sinhype_term1;
        denominator_x = 1.0/hd_x;

        exp_term_x = exp(-(gamma_u - l_2pi_kappa_sqrd)*t);
     
        sinhype_term2 =  sinh(l_2pi_over_b_m_sigma *(r-b))*r_r;

        f_x = exp_term_x * sinhype_term2 * denominator_x;

        
    }
    else{
        f_x = 0.0;
    }


    result = factor*(sum_y + f_x);


	return result;
}