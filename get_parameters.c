//
//  get_parameters.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/25.
//
//

#include "get_parameters.h"


void get_parameters_bimol1(std::vector<Species>& species_list, int I, int J, double* ka, double* D, double* sigma){

	*ka = 0.0;
	*D = 0.0;
	*sigma = 0.0;
	int i = 0;
	for (i = 0 ; i < species_list[I].asso_channel ; i++){
		if (sttoint(species_list[I].asso[i].partner.c_str() , species_list) == J){
			*ka = species_list[I].asso[i].rate_ka;
		}
	}
	*sigma = species_list[I].radius + species_list[J].radius;	
	
	*D = species_list[I].DiffCoeff + species_list[J].DiffCoeff;	
}


void get_parameters_bimol2(std::vector<Species>& species_list, int I, int J, double* ka, double* D, double* sigma, double* kd, double* gamma_b, double* gamma_u){

	*ka = 0.0; *D = 0.0; *sigma = 0.0; *kd = 0.0; *gamma_b = 0.0; *gamma_u = 0.0;

	*ka = species_list[I].disso[J].disso_rate.ka;

	*D = species_list[sttoint(species_list[I].disso[J].product1.c_str(), species_list)].DiffCoeff +species_list[sttoint(species_list[I].disso[J].product2.c_str(), species_list)].DiffCoeff;

	*sigma = species_list[sttoint(species_list[I].disso[J].product1.c_str(), species_list)].radius + species_list[sttoint(species_list[I].disso[J].product2.c_str(), species_list)].radius;

	*kd = species_list[I].disso[J].disso_rate.kd;

	*gamma_b = species_list[I].disso[J].disso_rate.gamma_b;

	*gamma_u = species_list[I].disso[J].disso_rate.gamma_u;
}

