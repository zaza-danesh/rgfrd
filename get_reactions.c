//	get_reactions.c

 
#include "get_reactions.h"

void get_reactions(int ns, std::vector<Species>& species_list, double ka1, double kd1, double k1, double ka2, double kd2, double k2){
	
	int i = 0;
	std::string letter;

	std::ifstream species_file;
	species_file.open(DATA_PATH SPECIES);
	for (i = 0 ; i < 6 ;i++){
		species_file >> letter;
	}

	Species new_species;
	new_species.name = "0";
	species_list.push_back(new_species);

	for (i = 1 ; i < ns ; i++){
		new_species.index = i;
		species_file >> new_species.name;
		species_file >> new_species.radius;
		species_file >> new_species.DiffCoeff;
		species_file >> new_species.quantity;
		species_list.push_back(new_species);
	}
	
species_file.close();


for (i = 1 ; i < ns ; i++){
	if (species_list[i].name == "S"){
	species_list[i].asso_flag = 1;
	species_list[i].asso_channel = 1;
	species_list[i].disso_flag = 0;
	species_list[i].disso_channel = 0;
	species_list[i].trans_monomol_flag = 0;
	species_list[i].trans_monomol_channel = 0;
	species_list[i].trans_bimol_flag = 0;
	species_list[i].trans_bimol_channel = 0;
	species_list[i].asso[0].partner = "E1";
	species_list[i].asso[0].product = "SE1";
	species_list[i].asso[0].rate_ka = ka1;
	species_list[i].total_transform_rate = 0.0;


}
else if (species_list[i].name == "P"){
	species_list[i].asso_flag = 1;
	species_list[i].asso_channel = 1;
	species_list[i].disso_flag = 0;
	species_list[i].disso_channel = 0;
	species_list[i].trans_monomol_flag = 0;
	species_list[i].trans_monomol_channel = 0;
	species_list[i].trans_bimol_flag = 0;
	species_list[i].trans_bimol_channel = 0;
	species_list[i].asso[0].partner = "E2";
	species_list[i].asso[0].product = "PE2";
	species_list[i].asso[0].rate_ka = ka2;
	species_list[i].total_transform_rate = 0.0;

}
else if (species_list[i].name == "E1"){
	species_list[i].asso_flag = 1;
	species_list[i].asso_channel = 1;
	species_list[i].disso_flag = 0;
	species_list[i].disso_channel = 0;
	species_list[i].trans_monomol_flag = 0;
	species_list[i].trans_monomol_channel = 0;
	species_list[i].trans_bimol_flag = 0;
	species_list[i].trans_bimol_channel = 0;
	species_list[i].asso[0].partner = "S";
	species_list[i].asso[0].product = "SE1";
	species_list[i].asso[0].rate_ka = ka1;
	species_list[i].total_transform_rate = 0.0;

}
else if (species_list[i].name == "E2"){
	species_list[i].asso_flag = 1;
	species_list[i].asso_channel = 1;
	species_list[i].disso_flag = 0;
	species_list[i].disso_channel = 0;
	species_list[i].trans_monomol_flag = 0;
	species_list[i].trans_monomol_channel = 0;
	species_list[i].trans_bimol_flag = 0;
	species_list[i].trans_bimol_channel = 0;
	species_list[i].asso[0].partner = "P";
	species_list[i].asso[0].product = "PE2";
	species_list[i].asso[0].rate_ka = ka2;
	species_list[i].total_transform_rate = 0.0;

}
else if (species_list[i].name == "SE1"){
	species_list[i].asso_flag = 0;
	species_list[i].asso_channel = 0;
	species_list[i].disso_flag = 1;
	species_list[i].disso_channel = 1;
	species_list[i].disso[0].product1 = "S";
	species_list[i].disso[0].product2 = "E1";
	species_list[i].disso[0].disso_rate.ka = ka1;
	species_list[i].disso[0].disso_rate.kd = kd1;
	species_list[i].disso[0].disso_rate.gamma_b = k1;
	species_list[i].disso[0].disso_rate.gamma_u = 0.0;
	species_list[i].trans_bimol_flag = 1;
	species_list[i].trans_bimol_channel = 1;
	species_list[i].trans_bimol[0].product1 = "P";
	species_list[i].trans_bimol[0].product2 = "E1";
	species_list[i].trans_bimol[0].trans_bimol_rate = k1;
	species_list[i].trans_bimol[1].product1 = "S";
	species_list[i].trans_bimol[1].product2 = "E1";
	species_list[i].trans_bimol[1].trans_bimol_rate = kd1;
	species_list[i].trans_monomol_flag = 0;
	species_list[i].trans_monomol_channel = 0;	
	species_list[i].total_transform_rate = k1 +kd1;
	
}
else if (species_list[i].name == "PE2"){
	species_list[i].asso_flag = 0;
	species_list[i].asso_channel = 0;
	species_list[i].disso_flag = 1;
	species_list[i].disso_channel = 1;
	species_list[i].disso[0].product1 = "P";
	species_list[i].disso[0].product2 = "E2";
	species_list[i].disso[0].disso_rate.ka = ka2;
	species_list[i].disso[0].disso_rate.kd = kd2;
	species_list[i].disso[0].disso_rate.gamma_b = k2;
	species_list[i].disso[0].disso_rate.gamma_u = 0.0;
	species_list[i].trans_bimol_flag = 1;
	species_list[i].trans_bimol_channel = 1;
	species_list[i].trans_bimol[0].product1 = "S";
	species_list[i].trans_bimol[0].product2 = "E2";
	species_list[i].trans_bimol[0].trans_bimol_rate = k2;
	species_list[i].trans_bimol[1].product1 = "P";
	species_list[i].trans_bimol[1].product2 = "E2";
	species_list[i].trans_bimol[1].trans_bimol_rate = kd2;
	species_list[i].trans_monomol_flag = 0;
	species_list[i].trans_monomol_channel = 0;	
	species_list[i].total_transform_rate = k2 + kd2;	
	}


}


}

