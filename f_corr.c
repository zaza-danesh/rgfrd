//
//  f_corr.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/06/13.
//
//
 

#include "f_corr.h"


double f_corr(double r, double t, phys_parameters_bimol1 physparams){

	double result = 0.0;
	const double ka = physparams.ka;
	const double D = physparams.D;
	const double sigma = physparams.sigma;
	const double r0 = physparams.r0;

	const double r_sqr = r*r;
	const double r_r0 = r*r0;
	const double sqrt_r_r0 = sqrt(r_r0);
	const double sqrt_r_r0_r  = 1.0/sqrt_r_r0;
	double	tolerance = 1e-8;

	gsl_integration_workspace*  workspace = gsl_integration_workspace_alloc(2000);

	result = -r_sqr*R0( r, t, workspace, tolerance, physparams);


	return result;

}




double R0 (double r, double t, gsl_integration_workspace* workspace, double tolerance, phys_parameters_bimol1 physparams ){

	double ka = physparams.ka;
	double D = physparams.D;
	double sigma = physparams.sigma;
	double r0 = physparams.r0;

	double integral = 0.0;
	double error = 0.0;
	double u_max = sqrt(40.0/D*t);
	tolerance = 1e-10;


    integrand_R0_parameters params = { r, t, ka, D, sigma, r0};

    gsl_function F;
    F.function = &integrate_R0;
    F.params = &params;

	gsl_integration_qag(&F, 0.0, u_max, tolerance, THETA_TOLERANCE, 2000, GSL_INTEG_GAUSS61, workspace, &integral, &error);


	return integral;
}

double integrate_R0(double u, void* params){

	struct integrand_R0_parameters *p  = (struct integrand_R0_parameters *) params;
	phys_parameters_bimol1 physparams;
	const double r = p->r;
	const double t = p->t;
	physparams.ka = p->ka;
	physparams.D = p->D;
	physparams.sigma = p->sigma;
	physparams.r0 = p->r0;
	
	
	return integrand_R0( u,  r, t, physparams);
}


double integrand_R0(double u, double r, double t, phys_parameters_bimol1 physparams){

	const double ka = physparams.ka;
	const double D = physparams.D;
	const double sigma = physparams.sigma;
	const double r0 = physparams.r0;

	const double real_n = 0.0;

	const double kD = 4*Pi*sigma*D;
	const double ka_sigma2 = 2*ka*sigma;
	const double ka_sigma2_over_kD = ka_sigma2/kD;
	const double ka_sigma2_over_kD_m_n2 = ka_sigma2_over_kD -2*real_n;

	const double u_sqr = u*u;
	const double term1 = exp(-D*t*u_sqr);

	const double s_u = sigma*u;
	const double r_u = r*u;
	const double r0_u = r0*u;
	// std::cout << "r" << r << std::endl;
	// const double js = spherical_j_bessel(0, s_u);

	const double js = gsl_sf_bessel_j0(s_u);
	// std::cout << js << std::endl;
	const double ys = spherical_y_bessel(0, s_u); 
	const double js1 = spherical_j_bessel(1, s_u);
	const double ys1 = spherical_y_bessel(1, s_u);
	const double jr = spherical_j_bessel(0, r_u);
	const double yr = spherical_y_bessel(0, r_u);
	const double jr0 = spherical_j_bessel(0, r0_u);
	const double yr0 = spherical_y_bessel(0, r0_u);


	const double R1 = (1 +ka_sigma2_over_kD_m_n2)*js + 2*s_u*js1;
	// std::cout << "R1 " << R1 << std::endl;
	const double R2 = (1 +ka_sigma2_over_kD_m_n2)*ys + 2*s_u*ys1;

	const double F1 = jr*jr0 - yr*yr0;
	const double F2 = jr0*yr + jr*yr0;

	const double F1R1 = R1*jr*jr0 - R1*yr*yr0;
	const double F2R2 = R2*jr0*yr + R2*jr*yr0;

	const double num = 2*u_sqr*R1*(F1R1 + F2R2);
	const double den = Pi*(R1*R1 + R2*R2);

	// const double result = term1*num/den;

	double result =u_sqr*exp(-u_sqr);

	return result;

}




