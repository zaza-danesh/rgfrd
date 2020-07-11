//
//  space_r_bimol2.c
//  
//
//  Created by Zahedeh Bashardanesh on 05/28/14.
//
//


 
#include "space_r_bimol2.h"
//This function gives the cumulative distribution function of space given the sampled time

double space_r_bimol2(double r,  double sampled_time, space_r_coeff_struct space_r_coeff){

    if ( r ==  space_r_coeff.b/ALPHA || sampled_time == 0.0){
      return 0.0;
    }

     int i = 0;
     int max_num_root = MAX_NUM_ROOT;


    double exp_term_x = 0.0;
    double exp_term_y = 0.0;
    double partial_sum = 0.0;
    double sum_y = 0.0;
    double f_x = 0.0;
    double result = 0.0;

    if (space_r_coeff.real_root_flag == 1){
        max_num_root--;
        f_x = space_r_coeff.denominator_x * (space_r_coeff.coshype_p_sinhype - space_r_coeff.coshype_sinhype_coeff * r *cosh(space_r_coeff.coshype_sinhype_arg*(r - space_r_coeff.b)) + space_r_coeff.coshype_sinhype_coeff * space_r_coeff.coshype_sinhype_coeff *  sinh(space_r_coeff.coshype_sinhype_arg *(r - space_r_coeff.b)));
        exp_term_x = exp(-space_r_coeff.exp_arg_x*sampled_time);
        f_x = exp_term_x * f_x;
    }
    else{
        f_x = 0.0;
    }



    for (i = 0; i < max_num_root ; i++){
        partial_sum = space_r_coeff.denominator_y[i] * (space_r_coeff.cosine_p_sine[i] - space_r_coeff.cosine_sine_coeff[i] * r * cos(space_r_coeff.cosine_sine_arg[i] * (r - space_r_coeff.b))  + space_r_coeff.cosine_sine_coeff[i] * space_r_coeff.cosine_sine_coeff[i] * sin(space_r_coeff.cosine_sine_arg[i] * (r - space_r_coeff.b)) );
        exp_term_y = exp(- space_r_coeff.exp_arg_y[i] * sampled_time);
        partial_sum = exp_term_y * partial_sum;
        sum_y = sum_y + partial_sum;
              
    }
    
    result = space_r_coeff.pre_factor*(sum_y - f_x);

    if ( fabs(result) > epsilon){
          return result;
     }
     else{
          return 0.0;
     }
    

}


