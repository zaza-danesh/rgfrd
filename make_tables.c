//
//  make_tables.c
//  
//
//  Created by Zahedeh Bashardanesh on 10/31/13.
//
//
 
#include "make_tables.h"

//**********************************************************************************************************************************
void makeRnTable(doubleVector& RnTable, double r, double t, phys_parameters_bimol1 physparams){

    RnTable.clear();
    unsigned int n = 0;
    const double ka = physparams.ka;
    const double D = physparams.D;
    const double sigma = physparams.sigma;
    const double r0 = physparams.r0;

     double absRn = 0.0;
     double Rn_value = 0.0;
     double Rn_prev = 0.0;
// std::cout << __LINE__ << std::endl;
    {
        const double p = p_irr(r, t, physparams);
        // std::cout << __LINE__ << std::endl;
        const double ip_free_max = ip_theta_free( Pi, r,t, physparams)*2*Pi*r*r;
// std::cout << __LINE__ << std::endl;
        if(fabs((p - ip_free_max)/ip_free_max) < 1e-8)
        {
            // std::cout << __LINE__ << std::endl;
            return;
        }
    }
// std::cout << __LINE__ << std::endl;
    const double pfree_max = p_free_max(r, t,  physparams);
    // std::cout << __LINE__ << std::endl;
    gsl_integration_workspace*  workspace = gsl_integration_workspace_alloc(2000);
    // std::cout << __LINE__ << std::endl;
    const double RnFactor = 1.0/(4.0*Pi* sqrt(r*r0));
    // std::cout << __LINE__ << std::endl;
    const double integrationTolerance = pfree_max / RnFactor * THETA_TOLERANCE;
    // std::cout << __LINE__ << std::endl;
    const double truncationTolerance = pfree_max * THETA_TOLERANCE * 1e-1;
    // std::cout << __LINE__ << std::endl;
    
// std::cout << __LINE__ << std::endl;
    for (;;) 
    {
// std::cout << __LINE__ << std::endl;
        Rn_value = Rn(n, r, t, workspace, integrationTolerance, physparams);
        // std::cout << __LINE__ << std::endl;
        RnTable.push_back(Rn_value);
        
        // truncate when converged enough.
        absRn = (fabs(Rn_value));
        if(absRn * RnFactor < truncationTolerance &&
            absRn < Rn_prev)
        {
            break;
        }

        if(n >= MAX_ORDER)
        {
            std::cout << "makeRnTable: Rn didn't converge" <<std::endl;
            break;
        }
        
        Rn_prev = fabs(Rn_value);
        
        ++n;
    }
    gsl_integration_workspace_free(workspace);
}


void makeLegendreTable(doubleVector const& RnTable, double theta ) {

   
    const int tableSize = RnTable.size();
    if(tableSize == 0)
    {
        std::cout << "tablesize = 0" << std::endl;

    }

    const double cos_theta = cos(theta);
    
    // lgndTable is offset by 1. lengTable[0] -> n = -1

    doubleVector lgndTable(tableSize + 2);
    lgndTable[0] = 1.0; // n = -1
    gsl_sf_legendre_Pl_array(tableSize, cos_theta, &lgndTable[1]);
}

//**********************************************************************************************************************************

double Rn(unsigned int n, double r, double t, gsl_integration_workspace* workspace, double tolerance, phys_parameters_bimol1 physparams ){

    double ka = physparams.ka;
    double D = physparams.D;
    double sigma = physparams.sigma;
    double r0 = physparams.r0;

    double integral = 0.0;
    double error = 0.0;
    double u_max = sqrt(40.0/D*t);


    integrand_Rn_parameters params = {n, r, t, ka, D, sigma, r0};

    gsl_function F;
    F.function = &integrate_Rn;
    F.params = &params;
    gsl_integration_qag(&F, 0.0, u_max, tolerance, THETA_TOLERANCE, 2000, GSL_INTEG_GAUSS61, workspace, &integral, &error);

 return integral;
}




//**********************************************************************************************************************************


double integrate_Rn(double u, void* params){

    struct integrand_Rn_parameters *p  = (struct integrand_Rn_parameters *) params;
    phys_parameters_bimol1 physparams;
    unsigned int n = p->n;
    double r = p->r;
    double t = p->t;
    physparams.ka = p->ka;
    physparams.D = p->D;
    physparams.sigma = p->sigma;
    physparams.r0 = p->r0;


    double R = integrand_Rn(n,  u,  r,  t, physparams);
    return R;
    
}


//**********************************************************************************************************************************


double integrand_Rn(unsigned int n, double u, double r, double t, phys_parameters_bimol1 physparams){

 const double ka = physparams.ka;
 const double D = physparams.D;
 const double sigma = physparams.sigma;
 const double r0 = physparams.r0;

 const double real_n = static_cast<double>(n);
 const double kD = 4*Pi*sigma*D;
 const double rr0 = r*r0;
 const double rr0_sqrd = sqrt(rr0);


 const double ka_over_kD = ka/kD; // = h is a dimensionless parameter you can find in the Jaeger and Carslaw
 const double ka_over_kD2 = 2.0*ka_over_kD;
 const double ka_over_kD2_m_2n = ka_over_kD2 - 2.0*real_n;

 const double u_sqr = u*u;
 const double term1 = exp(-D*t*u_sqr);

 const double s_u = sigma*u;
 const double r_u = r*u;
 const double r0_u = r0*u;

 const double js = spherical_j_bessel(n, s_u);
 const double ys = spherical_y_bessel(n, s_u); 
 const double js1 = spherical_j_bessel(n+1, s_u);
 const double ys1 = spherical_y_bessel(n+1, s_u);
 const double jr = spherical_j_bessel(n, r_u);
 const double yr = spherical_y_bessel(n, r_u);
 const double jr0 = spherical_j_bessel(n, r0_u);
 const double yr0 = spherical_y_bessel(n, r0_u);


 const double R1 = (1 + ka_over_kD2_m_2n)*js + 2*s_u*js1;
 const double R2 = (1 + ka_over_kD2_m_2n)*ys + 2*s_u*ys1;
 const double F1 = jr*jr0 - yr*yr0;
 const double F2 = jr0*yr + jr*yr0;
 const double F1R1 = R1*F1;
 const double F2R2 = R2*F2;

 const double R1_sqr = R1*R1;
 const double R2_sqr = R2*R2;

 const double num = 2*rr0_sqrd*u_sqr*R1*(F1R1 + F2R2);
 const double den = Pi*(R1_sqr + R2_sqr);

 const double result = term1*num/den; 
 return result;

}




//**********************************************************************************************************************************


double  p_irr(double r, double t, phys_parameters_bimol1 physparams){

    const double ka = physparams.ka;
    const double D = physparams.D;
    const double sigma = physparams.sigma;
    const double r0 = physparams.r0;

    const double sqrt_D = sqrt(D);
    const double kD = 4.0*Pi*sigma*D;
    const double alpha = (1.0 + (ka / kD)) * (sqrt_D/ sigma);

    const double Dt4 =4.0*D*t;
    const double r_p_r0_m_2sigma = r + r0 - 2.0 * sigma;

    const double num1 =exp(- gsl_pow_2(r - r0) / Dt4);
    const double num2 = exp(- gsl_pow_2(r_p_r0_m_2sigma) / Dt4);
    const double num3 = W(r_p_r0_m_2sigma / sqrt(Dt4), alpha * sqrt(t));

    const double num = (num1 + num2) / sqrt(4.0 * Pi * t) -  alpha * num3;
    const double den = (4.0 * Pi * r * r0 * sqrt_D);
    
    const double jacobian = (4.0 * Pi* r * r);
    
    double result = num / den;
    result = jacobian*result;
    

    return result;

}


//**********************************************************************************************************************************


double ip_theta_free(double theta, double r, double t, phys_parameters_bimol1 physparams){

 const double D = physparams.D;
 const double r0 = physparams.r0;

 const double Dt = D*t;
 const double Dt2 = Dt + Dt;
 const double rr0 = r*r0;
 const double rr0_over_2Dt = rr0/Dt2;
 const double rsq_p_r0sq_over_4Dt = (r*r + r0*r0)/(Dt2+Dt2);

 const double term1 = expm1(rr0_over_2Dt - rsq_p_r0sq_over_4Dt);
 const double term2 = expm1(rr0_over_2Dt*cos(theta) - rsq_p_r0sq_over_4Dt);
 const double den = 4.0*sqrt(Pi*Pi*Pi*Dt)*rr0;

 const double result = (term1-term2)/den;

 return result;
}


//**********************************************************************************************************************************

double p_free_max(double r , double t, phys_parameters_bimol1 physparams){

   //Z:  
    const double D = physparams.D;
    const double r0 = physparams.r0;
    const double Dt4 = 4*D*t;
    const double Dt4_r = 1.0/Dt4;
    const double sqr_Dt4_Pi = sqrt(Pi*Dt4);
    const double sqr_Dt4_Pi_r = 1.0/sqr_Dt4_Pi;

    const double term1 = exp(-gsl_pow_2(r-r0)*Dt4_r);
    const double term2 = sqr_Dt4_Pi_r*sqr_Dt4_Pi_r*sqr_Dt4_Pi_r;

    const double result = term1*term2;

    return result;
}
