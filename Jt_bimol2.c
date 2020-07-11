//
//  Jt_bimol2.c
//  
//
//  Created by Zahedeh Bashardanesh on 02/25/15.
//
//
 

#include "Jt_bimol2.h"


double Jt_bimol2(double t,  Jt_coeff_struct Jt_coeff){

    int i = 0;
    int max_num_root = MAX_NUM_ROOT;

    double exp_term = 0.0;
    double partial_sum = 0.0;
    double sum = 0.0;
    double result = 0.0;


    for (i = 0; i < max_num_root ; i++){
        exp_term = exp(- Jt_coeff.lambda1[i] * t);
        partial_sum = exp_term * Jt_coeff.lambda2[i] ;
        sum = sum + partial_sum;
    }

    result = Jt_coeff.pre_factor * sum;

    return result;
}




