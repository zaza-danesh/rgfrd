//
//  It_bimol2.c
//  
//
//  Created by Zahedeh Bashardanesh on 05/25/14.
//
//
 

#include "It_bimol2.h"


double It_bimol2(double t,  It_coeff_struct It_coeff ){

    // std::cout << "It_bimol2" << std::endl;
    if ( t == 0){
        return 0.0;
    }
    int i = 0;
    int max_num_root = MAX_NUM_ROOT;
    
    double exp_term = 0.0;
    double partial_sum = 0.0;
    double sum = 0.0;

    double result = 0.0;


    for (i = 0; i < max_num_root ; i++){
        exp_term = exp(-It_coeff.lambda1[i] * t);
        partial_sum = exp_term * It_coeff.lambda2[i];
        sum = sum + partial_sum;
        
    }

    result = It_coeff.pre_factor * sum;


    if (result > epsilon ) {
        return result;
    }
    else{
        return 0.0;
    }
}




