//
//  bound_or_unbound.c
//  
//
//  Created by Zahedeh Bashardanesh on 11/20/14.
//
//
 
#include "find_state.h"





std::string find_state_intermediate(double rnd, double sampled_time, phys_parameters_bimol2 physparams, double p1, double p2){

    
    if (!(rnd <1.0 && rnd > 0.0)) {
        std::cout << "Error: find_state::find_state_intermediate rnd is not a uniform distributed number in [0,1)" << std::endl;
    }
    
    std::string state = "--";

    double ka = physparams.ka;
    double D = physparams.D;
    double sigma = physparams.sigma;
    double b = physparams.b;
    double kd = physparams.kd;
    double gamma_b = physparams.gamma_b;
    double gamma_u = physparams.gamma_u;

     if (b < sigma){
        std::cout << "Error: drawTime::time_bound: b < sigma " << std::endl;
    }

    double probability = p1/(p1+p2);
    
    if (rnd < probability){
        state = "unbound"; //Particles are unbound
    }
    else {
        state = "bound"; //Particles are bound
    }

    return state;

}

std::string find_state_exit(double rnd, double sampled_time, phys_parameters_bimol2 physparams, double It, double Qt, double Jt){


    std::string state = "--";

    double ka = physparams.ka;
    double D = physparams.D;
    double sigma = physparams.sigma;
    double b = physparams.b;
    double kd = physparams.kd;
    double gamma_b = physparams.gamma_b;
    double gamma_u = physparams.gamma_u;
    double d_survival = 0.0;

    // double real_root = 0.0;
    // double char_roots[MAX_NUM_ROOT] = {0.0};


    if (physparams.b != physparams.sigma){
        d_survival = gamma_b*Qt + gamma_u*It + Jt;
    
        if (rnd < gamma_b*Qt/d_survival){
          state = "bound"; //transformaiton of bound state
        }
        else if (rnd < (gamma_u*It + gamma_b*Qt)/d_survival){
          state = "unbound"; //transformation of unbound state
        }
        else {
          state = "separation"; // separation of particles
        }
    }
    else{
        std::cout <<__FILE__ << SP << __LINE__ <<  ": Rrror in switch == 4 where b = sigma" << std::endl;
        state = "no_exit"; // b = sigma, this happens when the particle subject to dissociate is one of the particles in pair
    }    
return state;
}

