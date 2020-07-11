//
//  mystructure.h
//  
//
//  Created by Zahedeh Bashardanesh on 2016/01/21.
//
//

#ifndef MYSTRUCTURE_H
#define MYSTRUCTURE_H

#ifdef __linux__
    #include "/home/zahedeh/MyInclude/mystlib.h"
#else 
    #include "/Users/zahedeh/Dropbox/Sorin_Zahedeh/Implementation/MyInclude/mystlib.h"
#endif


// #define LIB_PATH "/mystlib.h"

#define DEBUG std::cout << std::endl; std::cout << __FILE__ << SP << __LINE__ << std::endl; std::cout << std::endl;
#define DEBUG_specific std::cout << __FILE__ << SP << __LINE__ << SP 
#define DEBUGK  if (k == 21){std::cout << std::endl; std::cout << __FILE__ << SP << __LINE__ << std::endl; std::cout << std::endl;}

#define ADJUSTMENT_STEPS_NUMBER 100
// #define epsilon_m 2e-4 
// #define epsilon_w 1e-4
#define P_D_F 1
#define ALPHA 16
#define MIN_DIFFUSION_TIME_STEP 1e-10
#define INNER_SHELL 1e-3
#define MIN_CUTOFF_TIME 1
#define MAX_SHORT_TIME_COUNTER 30



#define INF 1e18
#define INFTIME 1e100
#define MAX_NUM_ROOT 30
#define MAX_TIME_ALG 1e3 
// #define MAX_SHORT_TIME_COUNTER 30
#define epsilon 1e-8
#define rand_epsilon 1e-6
#define EPSILON  1e-8
#define tspacing .1 
#define FINAL_TIME_SKETCH 100
#define rspacing  .001
#define theta_spacing Pi/100

#define tmax 50
#define THETA_TOLERANCE  1e-5
#define MAX_ORDER 70
// #define PAIR_DIST 1e-2
#define INNER_SHELL 1e-3
#define MIN_TIME_STEP 1e-10

#define DATA_PATH ""
#define DATA_OUTPUT_PATH ""
#define SKETCH_PATH ""
#define SYSTEM 		"system.txt"
#define SPECIES 	"species.txt"
#define REACTIONS  "reactions.txt"

#define TIME_STEP_FILE_NAME "time_step.dat"
#define NUMBER_OF_PARTICLES_FILE_NAME "number_of_particles.dat"
#define B_VALUE_FILE_NAME "b.dat"
#define  PARTICLES_N_FILE_NAME "particles_n.dat"
#define  PARTICLES_N2_FILE_NAME "particles_n2.dat"
#define  PARTICLES_R_FILE_NAME "particles_r.dat"
#define  PARTICLES_R2_FILE_NAME "particles_r2.dat"

#define PARTICLES_X_FILE_NAME "particles_x.dat"
#define PARTICLES_Y_FILE_NAME "particles_y.dat"
#define PARTICLES_Z_FILE_NAME "particles_z.dat"

extern int KKKK ;

typedef double REAL;
struct TimeTable{

	double t;
	struct Particle;
};



struct Reactants {
	char type1;
	char type2;
};

struct Product {
	char type;
};

struct Reaction{
	char m_or_b;
	std::string type;
	std::string formula;
};

typedef std::vector<double> doubleVector;
typedef std::vector<TimeTable> timetable;
typedef std::vector<Reactants> ReactantsList;
typedef std::vector<Product> ProductList;


struct particle_particle_diffusion_time_struct{
	int pi_1;
	int pi_2;
	double t_m;
	double dist;
	particle_particle_diffusion_time_struct(){
		pi_1 = -1;
		pi_2 = -1;
		t_m = 0.0;
		dist = 0.0;
	}
};

struct particle_boundary_diffusion_time_struct{
	int pi;
	int si;
	double t_w;
	double dist;
	particle_boundary_diffusion_time_struct(){
		pi = -1;
		si = -1;
		t_w = 0.0;
		dist = 0.0;
	}
};




struct Partner{
	std::string partner;
	std::string product;
	double ka;
	double kd;
	double b;
	Partner(){
		partner = "-";
		product = "-";
		ka = 0.0;
		kd = 0.0;
		b = 0.0;
	}

};

struct Asso_Channel{
	std::string partner;
	std::string product;
	double rate_ka; //ka
	Asso_Channel(){
		partner = "-";
		product = "-";
		rate_ka = 0.0;
	}

};
struct Disso_Rate{
	double ka;
	double kd;
	double gamma_b;
	double gamma_u;
	Disso_Rate(){
		ka = 0.0;
		kd = 0.0;
		gamma_b = 0.0;
		gamma_u = 0.0;
	}
};

struct Disso_Channel{
	std::string product1;
	std::string product2;
	Disso_Rate disso_rate;
	Disso_Channel(){
		product1 = "-";
		product2 = "-";
	}

};

struct Trans_Bimol_Channel{
	std::string product1;
	std::string product2;
	double trans_bimol_rate;
	Trans_Bimol_Channel(){
		product1 = "-";
		product2 = "-";
		trans_bimol_rate = 0.0;
	}

};

struct Trans_Monomol_Channel{
	std::string product;
	double trans_monomol_rate;
	Trans_Monomol_Channel(){
		product = "-";
		trans_monomol_rate = 0.0;
	}

};


struct Species{
	std::string name;
	int index;
	int quantity;
	int asso_flag;
	int asso_channel;
	int disso_flag;
	int disso_channel;
	int trans_monomol_flag;
	int trans_monomol_channel;
	int trans_bimol_flag;
	int trans_bimol_channel;	
	double radius;
	double DiffCoeff;
	double total_transform_rate;
	Asso_Channel asso[5]; // the number of channels is limited to 3, it shouldn't be more than asso_channel
	Disso_Channel disso[5]; // the number of channels is limited to 3, it shouldn't be more than disso_channel
	Trans_Monomol_Channel trans_monomol[5]; // the number of channels is limited to 3, it shouldn't be more than trans_monomol_channel
	Trans_Bimol_Channel trans_bimol[5];  // the number of channels is limited to 3, it shouldn't be more than trans_bimol_channel
	


	Species(){
		name = "-";
		quantity = 0;
		asso_flag = 0;
		asso_channel = 0;
		disso_flag = 0;
		disso_channel = 0;
		trans_monomol_flag = 0;
		trans_monomol_channel = 0;
		trans_bimol_flag = 0;
		trans_bimol_channel = 0;
		radius = 0.0;
		DiffCoeff = 0.0;
		total_transform_rate = 0.0;

	}
};

struct Particle{
	int index;
	int pair_flag;
	int update_flag;
	int exist_flag;
	int reflect_flag;
	double t_w;
	double x;
	double y;
	double z;
	double r;
	double phi;
	double theta;
	double deltaX;
	double deltaY;
	double deltaZ;

	Particle(){
		index = -1;
		pair_flag = 0;		
		update_flag = 0;
		exist_flag = 1;
		reflect_flag = 0;
		t_w = 0.0;
		x = 0.0;
		y = 0.0;
		z = 0.0;
		r = 0.0;
		theta = 0.0;
		phi = 0.0;
		deltaX = 0.0;
		deltaY = 0.0;
		deltaZ = 0.0;		
	}
};

struct Roots{
	int index;
	double real_root;
	double char_roots[MAX_NUM_ROOT];
};


struct It_coeff_struct{
		int index;
		double pre_factor;
		double lambda1[MAX_NUM_ROOT];
		double lambda2[MAX_NUM_ROOT];
};

struct Qt_coeff_struct{
		int index;
		double pre_factor1;
		double pre_factor2;
		double lambda1[MAX_NUM_ROOT];
		double lambda2[MAX_NUM_ROOT];
};

struct Jt_coeff_struct{
		int index;
		double pre_factor;
		double lambda1[MAX_NUM_ROOT];
		double lambda2[MAX_NUM_ROOT];
};

struct space_r_coeff_struct{
	int index;
	int real_root_flag;
	double b;
	double pre_factor;
	double denominator_x;
	double coshype_p_sinhype;
	double coshype_sinhype_coeff;
	double coshype_sinhype_arg;
	double exp_arg_x;
	double denominator_y[MAX_NUM_ROOT];
	double cosine_p_sine[MAX_NUM_ROOT];
	double cosine_sine_coeff[MAX_NUM_ROOT];
	double cosine_sine_arg[MAX_NUM_ROOT];
	double exp_arg_y[MAX_NUM_ROOT];

};
struct pair_struct {
	int p1; 
	int p2;
	double dist;
	pair_struct(){
		p1 = -1;
		p2 = -1;
		dist = 0.0;
	} 
};

struct phys_parameters_bimol1{
	double ka;
	double D;
	double sigma;
	double r0;
	phys_parameters_bimol1(){
		ka = 0.0;
		D = 0.0;
		sigma = 0.0;
		r0 = 0.0;
	}
};

struct phys_parameters_bimol2{
	double ka;
	double D;
	double sigma;
	double b;
	double kd;
    double gamma_b;
    double gamma_u;
    phys_parameters_bimol2(){
    	ka = 0.0;
    	D = 0.0;
    	sigma = 0.0;
    	b = 0.0;
    	kd = 0.0;
    	gamma_b = 0.0;
    	gamma_u = 0.0;
    }
};

struct event_struct{
	std::string event;
	int pair_index;
	int pi_bimol2;
	int pi_monomol;
	double t_r;
	double physparams_monomol;
	struct phys_parameters_bimol1 physparams_bimol1;
	struct phys_parameters_bimol2 physparams_bimol2;
	event_struct(){
		event = "--";
		t_r = INF;
		pair_index = -1;
		pi_bimol2 = -1;
		pi_monomol = -1;
		physparams_monomol = 0.0;
		physparams_bimol1.ka = 0.0;
		physparams_bimol1.D = 0.0;
		physparams_bimol1.sigma = 0.0;
		physparams_bimol1.r0 = 0.0;
		physparams_bimol2.ka = 0.0; 
		physparams_bimol2.D = 0.0;
		physparams_bimol2.sigma = 0.0;
		physparams_bimol2.b = 0.0;
		physparams_bimol2.kd = 0.0;
		physparams_bimol2.gamma_b = 0.0;
		physparams_bimol2.gamma_u = 0.0;
	}
};

#endif
