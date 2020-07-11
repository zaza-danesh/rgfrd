//
//  f_free.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/06/13.
//
//
 

#include "f_free.h"





double f_free(double r, double t, phys_parameters_bimol1 physparams){

	double result = 0.0;
	const double ka = physparams.ka;
	const double D = physparams.D;
	// const double sigma = physparams.sigma;
	const double r0 = physparams.r0;
	const double r_sqr = r*r;
	const double r_p_r0 = r + r0;
	const double r_p_r0_sqr = r_p_r0*r_p_r0;
	const double r_r0 = r*r0;

	const double Dt = D*t;
	const double Dt4 = 4*Dt;
	const double Pi2 = 2*Pi;
	double den = sqrt(Pi2*Dt*Dt*Dt);
	den = den*r_r0;

	double factor = r_sqr/den;
	factor = 2*Dt*factor;

	const double num = exp(-r_p_r0_sqr/Dt4)*expm1(r_r0/Dt);
	

	result = factor*num;
	
	return result;
}


// double integrand_J0(double u, double r, double t, phys_parameters_association physparams){

// 	// const double ka = physparams.ka;
// 	const double D = physparams.D;
// 	const double sigma = physparams.sigma;
// 	const double r0 = physparams.r0;

	
// 	const double u_sqr = u*u;
// 	const double term1 = -exp(-D*t*u_sqr);

// 	const double s_u = sigma*u;
// 	const double r_u = r*u;
// 	const double r0_u = r0*u;

	
// 	const double jr = sj_bessel(0, r_u);
// 	const double jr0 = sj_bessel(0, r0_u);
	

// 	const double num = jr*jr0*u_sqr;

// 	const double result = term1*num;


// 	return result;

// }


// double integrate_J0(double u, void* params){

// 	struct integrand_J0_parameters *p  = (struct integrand_J0_parameters *) params;
// 	phys_parameters_association physparams;
// 	const double r = p->r;
// 	const double t = p->t;
// 	physparams.ka = p->ka;
// 	physparams.D = p->D;
// 	physparams.sigma = p->sigma;
// 	physparams.r0 = p->r0;

// 	return integrand_J0( u, r, t, physparams);
	
// }

// double J0( double r, double t, gsl_integration_workspace* workspace, double tolerance, phys_parameters_association physparams ){

// 	double ka = physparams.ka;
// 	double D = physparams.D;
// 	double sigma = physparams.sigma;
// 	double r0 = physparams.r0;

// 	double integral = 0.0;
// 	double error = 0.0;
// 	double u_max = sqrt(40.0/D*t);  // for J_n+0.5 I don't know what to set as the upper limit, 


//     integrand_J0_parameters params = {r, t, ka, D, sigma, r0};

//     gsl_function F;
//     F.function = &integrate_J0;
//     F.params = &params;

// 	gsl_integration_qag(&F, 0.0, umax, tolerance, THETA_TOLERANCE, 2000, GSL_INTEG_GAUSS61, workspace, &integral, &error);

// 	return integral;
// }