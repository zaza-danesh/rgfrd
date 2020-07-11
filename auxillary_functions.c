//
//  auxillary_functions.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/28/13.
//
//
 

#include "auxillary_functions.h"

double expsq_erfc(double x)
{
    double result;

    const double x_sq = x*x;
    const double x_r = 1.0/x;
    double term1 = 0.0;
	double term2 = 0.0;

    // Asymptotic expansion
    if(x > 26.0)
    {

    	const double Pi_sqrd = sqrt(Pi);
    	const double Pi_sqrd_r = 1.0/(Pi_sqrd);

        const double x2_sq_r = (1.0 / (2.0 * x_sq));  // 2 / (2 x)^2

    //       up to second term in the expansion.
    //       abs err ~= 9e-8 at x == 20, 3e-8 at x == 25

    //       the third term 
    //       - (8 / (x2sq * x2sq * x2sq))       
    //       and beyond doesn't have a major contribution for large x.
        

        term1 = (Pi_sqrd_r*x_r);
        term2 =  1 - x2_sq_r + x2_sq_r*x2_sq_r;

        result = term1*term2;
    }
     else
     {
        result = exp(x_sq) * erfc(x);
     }
    return result;
}

double  W(double a, double b ){

	  return (exp(- a * a) * expsq_erfc(a + b));

}

double spherical_j_bessel(unsigned int n, double x){

    if (n <= 3){
        return j_small_n(n,  x);
    }
    else{
     
     return gsl_sf_bessel_jl(n, x);    
    }
}

double spherical_y_bessel(unsigned int n, double x){
    
    if (n <= 3){
        // std::cout << y_small_n << std::endl;
        return y_small_n( n,  x);
    }
    else{
        return gsl_sf_bessel_yl(n, x);
    }
}

double j_small_n(unsigned int n, double x){

    if (n==0)
    {
        return sin(x)/x;
    }

    if (n==1)
    {
        double x_r = 1.0/x;
        double x_sq_r = x_r*x_r;
        return (sin(x)*x_sq_r  - cos(x)*x_r);
    }

    if( n==2)
    {
        double x_r = 1.0/x;
        double x_sq_r = x_r*x_r;
        double x_tr_r = x_sq_r*x_r;

        return ((3*x_tr_r - x_r)*sin(x) -3*x_sq_r*cos(x));
    }
    if( n==3)
    {
        double x_r = 1.0/x;
        double x_sq_r = x_r*x_r;
        double x_tr_r = x_sq_r*x_r;

        return ((15*x_sq_r*x_sq_r - 6*x_sq_r)*sin(x) -(15*x_tr_r - x_r)*cos(x));
    }   
    else{
        std::cout << " auxillary_functions: n not in range" << std::endl;
    }
}


double y_small_n(unsigned int n, double x){

    if (n==0)
    {
        
        return -cos(x)/x;
    }

    if (n==1)
    {
        double x_r = 1.0/x;
        double x_sq_r = x_r*x_r;
        return (-cos(x)*x_sq_r  - sin(x)*x_r);
    }
    if( n==2)
    {
        double x_r = 1.0/x;
        double x_sq_r = x_r*x_r;
        double x_tr_r = x_sq_r*x_r;

        return (-(3*x_tr_r - x_r)*cos(x) -3*x_sq_r*sin(x));
    }
    if( n==3)
    {
        double x_r = 1.0/x;
        double x_sq_r = x_r*x_r;
        double x_tr_r = x_sq_r*x_r;

        return ((-15*x_sq_r*x_sq_r + 6*x_sq_r)*cos(x) -(15*x_tr_r - x_r)*sin(x));
    }

    else{
        std::cout << " auxillary_functions : n not in range" << std::endl;
    }
}











