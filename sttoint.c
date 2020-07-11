//
//  sttoint.c
//  
//
//  Created by Zahedeh Bashardanesh on 03/14/14.
//
//
 
#include "sttoint.h"

int sttoint(const char name[], std::vector<Species>& species_list){	

	int index= 0;
		while (name != species_list[index].name){
			index++;
		}	
	return index;

}