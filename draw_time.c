//
//  draw_time.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/23/13.
//
//
 
#include "draw_time.h"

 

double draw_time_bimol1(double rnd,  phys_parameters_bimol1 physparams){

     
    if (!(rnd <1.0 && rnd > 0.0)) {
        std::cout << "Error: draw_time::time_bimol1: rnd is not a uniform distributed number in [0,1)" << std::endl;
    }


 
    const unsigned int maxIter = 120;
    unsigned int i = 0;
 
    double sampled_time = 0.0;
    double low = 1e-100;
    // double high =1e+5;
    double high = MAX_TIME_ALG;

    double ka = physparams.ka;
    double D = physparams.D;
    double sigma = physparams.sigma;
    double r0 = physparams.r0;
 
    double kD = 4*Pi*sigma*D;
    double ka_p_kD = ka + kD;
    double r0_ka_p_kD = r0*ka_p_kD;
    double r0_ka_p_kD_r = 1/r0_ka_p_kD;
    double factor = sigma*ka*r0_ka_p_kD_r;
  
  
      
    if (sigma > r0){
        std::cout << "Error: draw_time::time_bimol1: sigma > r0" << std::endl;
    }
  

    if( rnd > time_bimol1(MAX_TIME_ALG , ka, D, sigma, r0 )){
        return (MAX_TIME_ALG + rnd);
    }

    
    struct time_bimol1_parameters params = {ka, D, sigma, r0, rnd}; 
    gsl_function F;
    F.function = &time_bimol1_rnd;
    F.params = &params;

    const gsl_root_fsolver_type* solverType = gsl_root_fsolver_brent;
    gsl_root_fsolver* solver = gsl_root_fsolver_alloc(solverType);
    gsl_root_fsolver_set(solver, &F, low, high);
    for (;;){
        gsl_root_fsolver_iterate(solver);
        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-10, 1e-3));
        if(status == GSL_CONTINUE){
            if(i >= maxIter){
                gsl_root_fsolver_free(solver);
                /*throw runtime_error("draw_time: failed to converge");*/
                std::cout << "Error:  draw_time::time_bimol1 failed to converge" << std::endl;
            }
        }
        else{
            break;
        }
        ++i;
    }
    
    sampled_time = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);

      
    return sampled_time;
}



// //**********************************************************************************************************************************************

double draw_time_bimol2(int k, double rnd,  phys_parameters_bimol2 physparams, It_coeff_struct It_coeff, Qt_coeff_struct Qt_coeff){


    if (!(rnd <1.0 && rnd > 0.0)) {
        std::cout << "Error: draw_time::time_bimol2: rnd is not a uniform distributed number in [0,1): " << rnd <<   std::endl;
    }

    int ii=0;
    unsigned int maxIter = 110;    
    unsigned int i = 0;
    double sampled_time = 0.0;
    double low = 2e-100;
    double high = MAX_TIME_ALG;



    if (physparams.b < physparams.sigma){
        std::cout << "Error: draw_time::time_bimol2: b < sigma " << std::endl;
        std::exit;
    }

    if ( rnd > time_bimol2(MAX_TIME_ALG,  It_coeff, Qt_coeff) ){
        return (MAX_TIME_ALG + rnd);
    }

    struct time_bimol2_parameters params;
    params.It_coeff = It_coeff;
    params.Qt_coeff = Qt_coeff;
    params.rnd = rnd;

    gsl_function F;
    F.function = &time_bimol2_rnd;
    F.params = &params;

    const gsl_root_fsolver_type* solverType = gsl_root_fsolver_brent;
    gsl_root_fsolver* solver = gsl_root_fsolver_alloc(solverType);

    gsl_root_fsolver_set(solver, &F, low, high);
    

    for (;;){
        gsl_root_fsolver_iterate(solver);
        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-10, 1e-3));
        if(status == GSL_CONTINUE){
            if(i >= maxIter){
                gsl_root_fsolver_free(solver);
                /*throw runtime_error("draw_time: failed to converge");*/
                std::cout << "Error: draw_time::time_bimol2 failed to converge" << std::endl;
            }
        }
        else{
            break;
        }
        ++i;
    }      
    sampled_time = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);

    return sampled_time;
}


// //**********************************************************************************************************************************************

double draw_time_monomol(double rnd, double kd){

    if (!(rnd <1.0 && rnd > 0.0)){
       std::cout << "Error: draw_time::time_monomol: rnd is not a uniform distributed number in [0,1): " <<  rnd << std::endl;
    }

    const unsigned int maxIter = 100;
    unsigned int i = 0;
    double sampled_time = 0.0;
 
    double low = 1e-100;
    double high = MAX_TIME_ALG;

    if( rnd > time_monomol_distribution( MAX_TIME_ALG,  kd ) ){
        return (MAX_TIME_ALG + rnd);
    }

    struct time_monomol_params params = {kd, rnd}; 
    gsl_function F;
    F.function = &time_monomol_rnd;
    F.params = &params;
 
 
    const gsl_root_fsolver_type* solverType = gsl_root_fsolver_brent;
    gsl_root_fsolver* solver = gsl_root_fsolver_alloc(solverType);
    gsl_root_fsolver_set(solver, &F, low, high);
   
    for (;;){
        gsl_root_fsolver_iterate(solver);
 
        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-10, 1e-3));
 
        if(status == GSL_CONTINUE){
            if(i >= maxIter){
                gsl_root_fsolver_free(solver);
                /*throw runtime_error("draw_time: failed to converge");*/
                std::cout << "Error:  draw_time::time_bimol2 failed to converge" << std::endl;
            }
        }
        else{
            break;
        }
        ++i;
    }
 
    sampled_time = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);
 
 
    return sampled_time;
}

