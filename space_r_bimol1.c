//
//  space_r_bimol1.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/29/13.
//
//
 
#include "space_r_bimol1.h"

double space_r_bimol1(double r, double t , double ka, double D, double sigma, double r0){

      double result = 0.0;
    
     const double r0_r = 1.0/r0;
     const double Dt = D*t;
     const double kD = 4*Pi*sigma*D;
     const double kD_r = 1.0/kD;
     const double sqrt_D = sqrt(D);
     const double sigma_r = 1.0/sigma;
     const double ka_over_kD = ka*kD_r;
     const double alpha = (1.0 + ka_over_kD)*sqrt_D*sigma_r;
     const double ka_p_kD = ka + kD;
     const double ka_p_kD_r = 1.0/ka_p_kD;
     const double Dt4 = 4*Dt;
     const double sqrt_Dt4 = sqrt(Dt4);
     const double sqrt_Dt4_r = 1.0/sqrt_Dt4;
     const double ka_sigma2 = 2*ka*sigma;
     const double sqrt_t = sqrt(t);
     const double alpha_sqrt_t = alpha*sqrt_t;

     const double r_p_r0_m_2sigma_over_sqrt_Dt4 = (r+r0-2.0*sigma)*sqrt_Dt4_r;
     const double r_m_r0_over_sqrt_Dt4 = (r-r0)*sqrt_Dt4_r;
     const double r0_m_sigma_over_sqrt_Dt4 = (r0-sigma)*sqrt_Dt4_r;

     const double term1 = (sqrt(Dt/Pi))*(expm1(-gsl_pow_2(r_p_r0_m_2sigma_over_sqrt_Dt4)) - expm1(-gsl_pow_2(r_m_r0_over_sqrt_Dt4)));
     
 
     double term21 = erf(r_m_r0_over_sqrt_Dt4) + erf(r_p_r0_m_2sigma_over_sqrt_Dt4);
     term21 = r0*ka_p_kD*term21;

     double term22 = erf(r0_m_sigma_over_sqrt_Dt4) - erf(r_p_r0_m_2sigma_over_sqrt_Dt4);
     term22 = ka_sigma2*term22;

     const double term2 = term21 + term22;


     double term31 = W(r0_m_sigma_over_sqrt_Dt4, alpha_sqrt_t);
     term31 = ka*sigma*term31;

     double term32 = W(r_p_r0_m_2sigma_over_sqrt_Dt4, alpha_sqrt_t);
     term32 = (ka*r + kD*(r-sigma))*term32;

     const double term3 = term31 - term32;

     result = r0_r*(term1 + ka_p_kD_r*(0.5*term2 + term3));

     return result;
}
