
//  check_max_time_step.c
//  
//
//  Created by Zahedeh Bashardanesh on 05/11/15.
//
//
 
#include "check_max_time_step.h"

int check_max_time_step(double& max_time_step,  int& short_time_counter){

	int time_flag = 0;
		if (max_time_step > MAX_TIME_ALG){
            max_time_step = MAX_TIME_ALG;
        }
        if (max_time_step <  MIN_TIME_STEP){
            short_time_counter++;
        }
        else{
            short_time_counter = 0;
        }
        
        if(short_time_counter > MAX_SHORT_TIME_COUNTER ){
            std::cout << __FILE__ << SP << __LINE__ << ": short_time_counter > MAX_SHORT_TIME_COUNTER " <<std::endl;
            time_flag = 1;
        }
    return time_flag;
}