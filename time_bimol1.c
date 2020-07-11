//
//  time_bimol1.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/28/13.
//
//

//Z: the time to the next reaction; survival probability -1 : S(t)-1
 
#include "time_bimol1.h"


double time_bimol1(double t, double ka, double D, double sigma,  double r0 ){
     double result = 0.0;
     const double sqrt_t = sqrt(t);
     const double sqrt_D = sqrt(D);
     const double sqrt_Dt = sqrt_t*sqrt_D;
     const double sqrt_Dt_r = 1.0/sqrt_Dt;

     const double kD = 4*Pi*sigma*D;
     const double kD_r = 1.0/kD;
     const double sigma_r = 1.0/sigma;
     const double ka_over_kD = ka*kD_r;
     const double alpha = (1.0 + ka_over_kD)*sqrt_D*sigma_r;
     const double ka_p_kD = ka + kD;
     const double r0_ka_p_kD = r0*ka_p_kD; 
     const double r0_ka_p_kD_r = 1/r0_ka_p_kD;
    
    const double arg1 = (r0-sigma)*0.5*sqrt_Dt_r;

    const double arg2 = alpha*sqrt_t;
    const double term1 = erfc(arg1);

    const double term2 = W(arg1, arg2);
     
    const double factor = sigma*ka*r0_ka_p_kD_r;


    result = factor*(term1 - term2);

    return result;
}



//**********************************************************************************************************************************************

double binding_probability(double t, double ka, double D, double sigma,  double r0 ){

    struct phys_parameters_bimol1 physparams;
    physparams.ka = ka;
    physparams.D = D;
    physparams.sigma = sigma;
    physparams.r0 = r0;
    double result = 0.0;


    const double sqrt_D = sqrt(D);
    const double sqrt_t = sqrt(t);
    const double sigma_r = 1.0/sigma;
    const double kD = 4*Pi*sigma*D;
    const double kD_r = 1.0/kD;
    const double ka_over_kD = ka*kD_r;
    const double ka_p_kD = ka + kD;
    const double r0_ka_p_kD = r0*ka_p_kD; 
    const double r0_ka_p_kD_r = 1/r0_ka_p_kD;
    const double factor = sigma*ka*r0_ka_p_kD_r; 


    const double arg1 = (r0-sigma)/(2*sqrt_D);
    const double arg2 = (1.0 + ka_over_kD)*sqrt_D*sigma_r;

    const double Term1 = exp(-4*arg1*arg2 - 3*arg1*arg1/t - 2*arg2*arg2*t);
    const double Term2 = sqrt_t*(arg1*(expm1(2*(arg1+arg2*t)*(arg1+arg2*t)/t)) + arg2*t);
    const double Term3 = sqrt(Pi)*(2*arg1*arg1 - arg2*arg2*t*t)*exp((arg1+arg2*t)*(arg1+arg2*t)/t)*erfc((arg1+arg2*t)/sqrt_t);
    const double num = Term1*(Term2-Term3);
    const double den = sqrt(Pi)*t*t;
    result = factor*num/den;
    return result;
   
}























