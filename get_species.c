//
//  get_species.c
//  
//
//  Created by Zahedeh Bashardanesh on 01/17/14.
//
//
 
#include "get_species.h"


int get_number_of_species(){
	std::string letter;
	int number_of_species = -1;
    std::ifstream species_file;
    
    species_file.open(DATA_PATH SPECIES);
    species_file >> letter;
    species_file >> number_of_species;
    number_of_species++;
    species_file.close();

    return number_of_species;
}





void get_species(int number_of_species, Species* species_list){

	int i;
	std::string rdnt;


	//open the file in read mode.
	std::ifstream infile;
	infile.open(DATA_PATH SPECIES);

	for (i = 0; i<7 ; i++){
		infile >> rdnt;
	}

	std::cout << std::endl;
	
	// std::cout <<  "Name:" << "		" <<"radius:" <<"		" <<  "Diff_Coeff:"<<"	"<<  "quantity:" <<std::endl;
	// 	std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
	for (i = 0 ; i< number_of_species ; i++){
		species_list[i].index = i;
		// std::cout << species_list[i].index << std::endl;
		infile >> species_list[i].name;
		std::cout << "name: " << species_list[i].name << SP;
		infile >> species_list[i].radius;			
		std::cout << "radius: " << species_list[i].radius << SP;
		infile >> species_list[i].DiffCoeff;
		std::cout << "D: " << species_list[i].DiffCoeff << SP;
		infile >> species_list[i].quantity;
		std::cout << "quantity: " << species_list[i].quantity << std::endl;


	}



	// std::cout << std::endl;
	// std::cout << "*********************************************************************************************" << std::endl;
	// std::cout << std::endl;

	infile.close();

}