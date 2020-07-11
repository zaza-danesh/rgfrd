//
//  draw_space.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/30/13.
//
//
 

#include "draw_space.h"



double draw_space_r_bimol1(double rnd, double t, phys_parameters_bimol1 physparams){

 if (!(rnd <1.0 && rnd > 0.0)) {
        std::cout << "draw_space::draw_space_r_bimol1: rnd is not a uniform distributed number in [0,1)" << std::endl;
    }



    const unsigned int maxIter = 100;
    unsigned int i = 0;
    unsigned int H = 0;

    double sampled_r = 0.0;

    double ka = physparams.ka;
    double D = physparams.D;
    double sigma = physparams.sigma;
    double r0 = physparams.r0;

    
    double kD = 4*Pi*sigma*D;
    double kD_r = 1.0/kD;
    double sqrt_D = sqrt(D);
	double sigma_r = 1.0/sigma;
	double ka_over_kD = ka*kD_r;
	double alpha = (1.0 + ka_over_kD)*sqrt_D*sigma_r;

	double low = r0; // should be adjusted from r0 to avoid root finding in the long tail  ********
    double high = r0;     // should be adjusted from r0 to avoid root finding in the long tail  ********
    const double sqrt_6Dt = sqrt(6*D*t);
    double value = 0.0;


    
   if (sigma > r0){
        std::cout << "draw_space::draw_space_r_bimol1: sigma > r0" << std::endl;
	}

	if (!(t>=0.0) ){
        std::cout << "draw_space::draw_space_r_bimol1: sampled_time negative! " << std::endl;
    } 

    if(t==0.0){
    	return r0;
    }

const double psurv = survival_bimol1(t, ka, D, sigma, r0);
rnd = rnd*psurv;
struct space_r_bimol1_parameters params = {ka, D, sigma, r0, rnd, t};



 	gsl_function F;
    F.function = &space_r_bimol1_rnd;
    F.params = &params;

// *********      *********      *********      *********      *********      *********      *********      *********      *********      
// const double sqrt6Dt(sqrt(6.0 * D * t));
    if(GSL_FN_EVAL(&F, r0) < 0.0)
    {
        // low = r0
        H = 3;

        for (;;)
        {
            
            high = r0 + H * sqrt_6Dt;
            // std::cout << "high: " << high << std::endl;
            value =GSL_FN_EVAL(&F, high);
            // std::cout << "value: " << value << std::endl;

            if(value > 0.0)
            {
                break;
            }

            ++H;

            if(H > 20)
            {
                // throw runtime_error("drawR: H > 20 while adjusting upper bound of r");
                std::cout << "draw_space::draw_space_r_bimol1: H>20 while adjusting upper bound of r" << std::endl;
                break;
            }
        }

    }
    else // GSL_FN_EVAL(&F, r0) > 0.0
    {
        // high = r0
         H = 3;

        for (;;)
        {
            low = r0 - H * sqrt_6Dt;
            if(low < sigma)
            {
                if(GSL_FN_EVAL(&F, sigma) > 0.0)
                {
                    // log_.info("drawR: p_int_r(sigma) > 0.0. "
                    //           "returning sigma.");
                    std::cout << "draw_space::draw_space_r_bimol1: displacement_irr(sigma)>0.0, returning sigma" << std::endl;
                    return sigma;
                }

                low = sigma;
                break;
            }

             value = GSL_FN_EVAL(&F, low);
            if(value < 0.0)
            {
                break;
            }

            ++H;
            if(H > 20)
            {
                // throw runtime_error("drawR: H > 20 while adjusting lower bound of r");
                std::cout << "draw_space::draw_space_r_bimol1: H>20 while adjusting lower bound of r" << std::endl;
                break;
            }
        }
    }

// *********      *********      *********      *********      *********      *********      *********      *********      *********      



const gsl_root_fsolver_type* solverType = gsl_root_fsolver_brent;

gsl_root_fsolver* solver = gsl_root_fsolver_alloc(solverType);

gsl_root_fsolver_set(solver, &F, low, high);


for (;;)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-5, 1e-2));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                /*throw runtime_error("drawDisplacement: failed to converge");*/
                std::cout << "draw_space::draw_space_r_bimol1: failed to converge" << std::endl;
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
     sampled_r = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);

return sampled_r;
}



// ***********************************************************************************************************************************************************************************  

double draw_space_theta_bimol1(double rnd, double r,  double t, phys_parameters_bimol1 physparams){
        if (!(rnd <1.0 && rnd > 0.0)) {
            std::cout << "drawTime: rnd is not a uniform distributed number in [0,1)" << std::endl;
        }        
        
    const unsigned int maxIter = 100;
    unsigned int i = 0;
    unsigned int H = 0;
    double low = 0;
    double high =Pi;

    double sampled_theta = 0.0;

    double ka = physparams.ka;
    double D = physparams.D;
    double sigma = physparams.sigma;
    double r0 = physparams.r0;

    // std::cout << ka << SP << D << SP << sigma << SP << r0 << SP << t << SP << r << std::endl;
    // std::cout << "ka: " << ka << SP << "D: " << D << SP << "sigma: " << sigma << SP << "r0: " << r0 << std::endl;

    
    double kD = 4*Pi*sigma*D;
    double kD_r = 1.0/kD;
    double sqrt_D = sqrt(D);
    double sigma_r = 1.0/sigma;
    double ka_over_kD = ka*kD_r;
    double alpha = (1.0 + ka_over_kD)*sqrt_D*sigma_r;
    if (!(r>=sigma))
    {
        std::cout << "drawTheta: sampled_r < sigma " << std::endl;    
    }

    if (!(r0>=sigma))
    {
        std::cout << "drawTheta: r0 < sigma " << std::endl;   
    }

    if (!(t>=0.0) ){
        std::cout << "drawTheta: sampled_timed negative! " << std::endl;
    } 


    if(t==0.0)
    { //no move

        return 0.0;
    }

    const double p_free_value = pp_free(Pi, r,t, physparams);
    const double p_corr_value = pp_corr(Pi, r,t, physparams);
    const double factor = p_free_value + p_corr_value;
    

    rnd = factor *rnd;

    struct space_theta_bimol1_parameters params = {ka, D, sigma, r0, rnd, t , r};
    gsl_function F;
    F.function = &space_theta_bimol1_rnd;    
    F.params = &params;
 
    const gsl_root_fsolver_type* solverType = gsl_root_fsolver_brent;

    gsl_root_fsolver* solver = gsl_root_fsolver_alloc(solverType);

    gsl_root_fsolver_set(solver, &F, low, high);

   
    for (;;)
    {

        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-5, 1e-2));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                /*throw runtime_error("drawTime: failed to converge");*/
                std::cout << " draw_space::space_theta_bimol1_rnd failed to converge" << std::endl;
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
    sampled_theta = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);

    return sampled_theta;
}


// ***********************************************************************************************************************************************************************************


double draw_space_r_bimol2(double rnd, double t,  phys_parameters_bimol2 physparams, space_r_coeff_struct space_r_coeff, double psurv){

    // std::cout << __FILE__ << std::endl;

    if (!(rnd <1.0 && rnd > 0.0)) {
        std::cout << "draw_space::draw_space_r_bimol2: rnd is not a uniform distributed number in [0,1)" << std::endl;
    }
    int ii=0;
    const unsigned int maxIter = 100;
    unsigned int i = 0;
    double sampled_r = 0.0;

       

    double ka = physparams.ka;
    double D = physparams.D;
    double sigma = physparams.sigma;
    double b = physparams.b;
    double kd = physparams.kd;
    double gamma_u = physparams.gamma_u;
    double gamma_b = physparams.gamma_b;

    double low = sigma ;
    double high = b ;


    // double real_root = 0.0;
    // double char_roots[MAX_NUM_ROOT] = {0.0};


    // find_roots(&real_root, char_roots,  physparams);
    
   if (b < sigma){
        std::cout << "draw_space::draw_space_r_bimol2: b < sigma " << std::endl;
    }

    if (!(t>=0.0) ){
        std::cout << "draw_space::draw_space_r_bimol2: sampled_time negative! " << std::endl;
    } 

if(t==0.0){

        return 0.0;
    }

// const double psurv = It_bimol2(t, ka, D, sigma, b, kd, gamma_b, gamma_u, species_roots.real_root, species_roots.char_roots);
rnd = rnd*psurv;
struct space_bimol2_parameters params;
params.space_r_coeff = space_r_coeff;
params.sampled_time = t;
params.rnd = rnd;


    gsl_function F;
    F.function = &space_r_bimol2_rnd;
    F.params = &params;

const gsl_root_fsolver_type* solverType = gsl_root_fsolver_brent;

gsl_root_fsolver* solver = gsl_root_fsolver_alloc(solverType);

gsl_root_fsolver_set(solver, &F, low, high);

for (;;)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-5, 1e-2));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                /*throw runtime_error("drawDisplacement: failed to converge");*/
                std::cout << "draw_space::draw_space_r_bimol2: failed to converge" << std::endl;
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
     sampled_r = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);

// std::cout << "Exit draw_space: bimol2 " << std::endl;
return sampled_r;

}









