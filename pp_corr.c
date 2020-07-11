//
//  pp_corr.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/06/13.
//
//
 

#include "pp_corr.h"

double pp_corr(double theta, double r, double t, phys_parameters_bimol1 physparams){

	int i = 0;
	double result = 0.0;
	const double ka = physparams.ka;
	const double D = physparams.D;
	const double sigma = physparams.sigma;
	const double r0 = physparams.r0;

	const double r_r0 = r*r0;
	const double sqr_r_r0 = sqrt(r_r0);
	const double Pi4 = 4*Pi;
	const double sqr_r_r0_Pi4 = Pi4*sqr_r_r0;
	const double factor = -1.0/sqr_r_r0_Pi4;
	const double cos_theta = cos(theta);

	double lgnd_n_m1 = 0.0;	
	double lgnd_n_p1 = 0.0;
	double Rn_Legendre = 0.0;
	double sum  = 0.0;

	doubleVector RnTable;

	for (i=0 ; i <RnTable.size() ; i++){
		RnTable[i] = 0.0;
	}

	makeRnTable(RnTable, r, t, physparams);
	const int  tableSize = RnTable.size();
	
    if(tableSize == 0)
    {
        return 0.0;
    }
    doubleVector lgndTable(tableSize + 2);

    lgndTable[0] = 1.0; // n = -1
    gsl_sf_legendre_Pl_array(tableSize, cos_theta, &lgndTable[1]);    

 	for (i = 0 ; i < tableSize ; i++){
 		lgnd_n_m1 = lgndTable[i];
 		lgnd_n_p1 = lgndTable[i+2];
 		Rn_Legendre = RnTable[i]*(lgnd_n_m1 - lgnd_n_p1);
 		
 		sum = sum + Rn_Legendre;
 	}
 	result = factor*sum;
	return result; 
}






double pp_corr_test(int k, double theta, double r, double t, phys_parameters_bimol1 physparams){

	int i = 0;
	double result = 0.0;
	const double ka = physparams.ka;
	const double D = physparams.D;
	const double sigma = physparams.sigma;
	const double r0 = physparams.r0;

	const double r_r0 = r*r0;
	const double sqr_r_r0 = sqrt(r_r0);
	const double Pi4 = 4*Pi;
	const double sqr_r_r0_Pi4 = Pi4*sqr_r_r0;
	const double factor = -1.0/sqr_r_r0_Pi4;
	const double cos_theta = cos(theta);

	double lgnd_n_m1 = 0.0;	
	double lgnd_n_p1 = 0.0;
	double Rn_Legendre = 0.0;
	double sum  = 0.0;

	doubleVector RnTable;

	for (i=0 ; i <RnTable.size() ; i++){
		RnTable[i] = 0.0;
	}
	makeRnTable(RnTable, r, t, physparams);
	const int  tableSize = RnTable.size();
	
    if(tableSize == 0)
    {
        return 0.0;
    }
    doubleVector lgndTable(tableSize + 2);
    lgndTable[0] = 1.0; // n = -1
    gsl_sf_legendre_Pl_array(tableSize, cos_theta, &lgndTable[1]);    
 	for (i = 0 ; i < tableSize ; i++){
 		lgnd_n_m1 = lgndTable[i];
 		lgnd_n_p1 = lgndTable[i+2];
 		Rn_Legendre = RnTable[i]*(lgnd_n_m1 - lgnd_n_p1);
 		
 		sum = sum + Rn_Legendre;
 	}
 	result = factor*sum;
	return result; 
}





