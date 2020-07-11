//
//  Qt_bimol2.c
//  
//
//  Created by Zahedeh Bashardanesh on 05/26/14.
//
//
 

#include "Qt_bimol2.h"


double Qt_bimol2(double t,  Qt_coeff_struct Qt_coeff){
    
    if (t == 0){
        return 1.0;
    }

    int i = 0;
    int max_num_root = MAX_NUM_ROOT;
    double expm1_term = 0.0;
    double partial_sum = 0.0;
    double sum = 0.0;
    double result = 0.0;

    
    for (i = 0; i < max_num_root ; i++){        
        expm1_term = expm1(-Qt_coeff.lambda1[i]*t);
        partial_sum = expm1_term * Qt_coeff.lambda2[i];
        sum = sum + partial_sum;              
    }

    result =  exp(-Qt_coeff.pre_factor2*t) * (1.0 + Qt_coeff.pre_factor1 * sum);

    return result;
}




