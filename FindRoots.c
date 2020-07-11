
//
//  Find_Roots.c
//  
//
//  Created by Zahedeh Bashardanesh on 2016/03/11.
//
//

#include "FindRoots.h"


std::vector<Roots> Find_Roots(std::vector<Species>& species_list, double b_value){


	double real_root = 0.0;
    double char_roots[MAX_NUM_ROOT] = {0.0};
	double ka = 0.0; double D = 0.0; double sigma = 0.0; double b = 0.0; double kd = 0.0; double gamma_b = 0.0; double gamma_u = 0.0;
	int J = 0;
	int i = 0;
	struct phys_parameters_bimol2 physparams;
	struct Roots roots;
	std::vector<Roots> species_roots;
	for (i = 0 ; i < species_list.size() ; i++){
		if (species_list[i].disso_flag !=0){
			roots.index = i;
			ka = species_list[i].disso[J].disso_rate.ka;
			D = species_list[sttoint(species_list[i].disso[J].product1.c_str(), species_list)].DiffCoeff +species_list[sttoint(species_list[i].disso[J].product2.c_str(), species_list)].DiffCoeff;
			sigma = species_list[sttoint(species_list[i].disso[J].product1.c_str(), species_list)].radius + species_list[sttoint(species_list[i].disso[J].product2.c_str(), species_list)].radius;
			// b = std::min( ALPHA*sigma, b_value + sigma);
			b = ALPHA*sigma;
			kd = species_list[i].disso[J].disso_rate.kd;
			gamma_b = species_list[i].disso[J].disso_rate.gamma_b;
			gamma_u = species_list[i].disso[J].disso_rate.gamma_u;
			physparams = make_parameters_bimol2(ka, D, sigma, b , kd, gamma_b, gamma_u );
			find_roots(roots.real_root, roots.char_roots,  physparams);
			species_roots.push_back(roots);
		}

	}
	return species_roots;
}


	