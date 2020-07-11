//
//  find_roots.c
//  
//
//  Created by Zahedeh Bashardanesh on 12/13/13
//  Modified by Zahedeh Bashardanesh on 11/04/14
//
 
#include "find_roots.h"


void find_roots( double& real_root, double* char_roots , phys_parameters_bimol2 physparams){

int i=0;
int j=0;
int n=0;
int k=0;
int h=0;
int m=0;
int real_root_flag = 0;
int multiple_root = 0;
int max_num_root = MAX_NUM_ROOT;

const double ka = physparams.ka;
const double D = physparams.D;
const double sigma = physparams.sigma;
const double b = physparams.b;
const double kd = physparams.kd; 
const double gamma_b = physparams.gamma_b;
const double gamma_u = physparams.gamma_u;
std::cout.precision(10);
const double b_m_sigma = b-sigma;
const double b_m_sigma_r = 1.0/b_m_sigma;
const double sqrt_D = sqrt(D);
const double sqrt_D_r = 1.0/sqrt_D;
const double kappa1_r = sigma*sqrt_D_r; 
const double kappa2_r = b*sqrt_D_r;
const double kappa_r = kappa2_r - kappa1_r;
const double kappa1 = 1./kappa1_r;
const double kappa2 = 1./kappa2_r;
const double kappa = 1./kappa_r;
const double kD = 4*Pi*sigma*D;
const double kD_r = 1.0/kD;
const double N_disso = 1 + ka*kD_r;
const double decay_diff = gamma_u-gamma_b;
double x_low = 0.0;
double x_high = 0.0;
double q_low = 0.0;
double q_high = 0.0;
double Delta = 0.0;
double value = 0.0;
double q = EPSILON;
double q_spacing = .01;
// double hy[1000] = {0};

    // std::ofstream hy_file;
    // hy_file.open(SKETCH_PATH "hy.dat");
    // for (i = 0 ; i<1000 ; i++){
    //   hy[i] = sketch_hy(q, kappa, kappa1_r, kd, N_disso, gamma_b, gamma_u);
    //   q = q + q_spacing;
    //   hy_file << hy[i] << std::endl;
    // }
    // hy_file.close();


if (N_disso > 1){
    double A = decay_diff-kd;
    double B = kappa1*N_disso;
    double C = B*decay_diff - kappa1*kd;
    double asymptote = 0.0;
    double fp_zero = 0.0;
    asymptote = kappa_r*sqrt(-C/B);
    asymptote = asymptote/(2*Pi);
    fp_zero = kappa*A/C; 
    Delta = N_disso*decay_diff - kd + kappa*kappa1_r*(decay_diff-kd);


    if (Delta > 0){
        max_num_root--;
        real_root = find_roots_hx(physparams, EPSILON, sqrt(gamma_u));
        if (fp_zero < 0 && fabs(fp_zero) >1){
            std::cout << "Error: Delta>0 and fp_zero < 0 && fabs(fp_zero) >1!"<<std::endl;
        }
        for (n = 1 ; n<max_num_root+1 ; n++){
            q_low = (2*n-1)/4.0;
            q_high = (2*n+1)/4.0;
            char_roots[n-1] = find_roots_hy(physparams, q_low, q_high);
        }
    }
    else { //Delta<0 
      real_root = 0.0;
        if (C<0){ //there is an asymptote
            if (asymptote < 0.25){//this works
                n = 0;
                char_roots[n] = find_roots_hy(physparams, asymptote, 0.25);
                for (n = 1 ; n<max_num_root ; n++){ //this works
                    q_low = (2*n-1)/4.0;
                    q_high = (2*n+1)/4.0;
                    char_roots[n] = find_roots_hy(physparams, q_low, q_high);
                }
            }
            else{ //if asymptote >0.25
                n = 1;
                while((2*n+1)/4.0 < asymptote){ //This works
                    q_low = (2*n-1)/4.0;
                    q_high = (2*n+1)/4.0;
                    char_roots[n-1] = find_roots_hy(physparams, q_low, q_high);
                    n++;
                }
                q_low = (2*n-1)/4.0;
                q_high = asymptote;
                char_roots[n-1] = find_roots_hy(physparams, q_low, q_high);
                q_low = asymptote;
                q_high = (2*n+1)/4.0;
                char_roots[n] = find_roots_hy(physparams, q_low, q_high);
                for(k = n+1; k<max_num_root ; k++){
                    q_low = (2*k-1)/4.0;
                    q_high = (2*k+1)/4.0;
                    char_roots[k] = find_roots_hy(physparams, q_low, q_high);
                }
            }
        }
        else{ //there is no asymptote //This works!
            n = 0;
            char_roots[n] = find_roots_hy(physparams, EPSILON, 0.25);
            for(n=1 ; n<max_num_root ; n++){
                q_low = (2*n-1)/4.0;
                q_high = (2*n+1)/4.0;
                char_roots[n] = find_roots_hy(physparams, q_low, q_high);
            }
        }
    }
}

else{ //N_disso = 1
        double arg = kd - (gamma_u - gamma_b);
    // std::cout << __LINE__ << "arg: " << arg << std::endl;
    if (arg>0){
        char_roots[0] = sqrt(arg)*kappa_r/(2*Pi);
        // std::cout << "char_roots[0]: " << char_roots[0] << std::endl;
        for (n = 1 ; n < max_num_root ; n++){
            q_low = (2*n-1)/4.0;
            q_high = (2*n+1)/4.0;
            char_roots[n] = find_roots_htan(physparams, q_low, q_high);
        }
    }
    else{
        for (n = 0 ; n < max_num_root ; n++){
            // std::cout << "n: " << n << std::endl;
            q_low = (2*n+1)/4.0;
            q_high = (2*n+3)/4.0;
            char_roots[n] = find_roots_htan(physparams, q_low, q_high);
        }    

    }
}


/*std::ofstream imaginary_roots_file;
imaginary_roots_file.open(SKETCH_PATH "imaginary_roots.dat");
imaginary_roots_file.precision(10);
for (i=0 ; i <max_num_root ; i++){
  imaginary_roots_file << char_roots[i] << std::endl;
}
// std::cout << __FILE__ <<SP << __LINE__ << SP << "root 0: " << char_roots[0] << std::endl;
imaginary_roots_file.close();*/

return;


}

double find_roots_hx(phys_parameters_bimol2 physparams, double x_low, double x_high){


const double ka = physparams.ka;
const double D = physparams.D;
const double sigma = physparams.sigma;
const double b = physparams.b;
const double kd = physparams.kd;
const double gamma_b = physparams.gamma_b;
const double gamma_u = physparams.gamma_u;

const double sqrt_D = sqrt(D);
const double sqrt_D_r = 1.0/sqrt_D;
const double kappa1_r = sigma*sqrt_D_r; 
const double kappa2_r = b*sqrt_D_r;
const double kappa_r = kappa2_r - kappa1_r;
const double kappa1 = 1./kappa1_r;
const double kappa2 = 1./kappa2_r;
const double kappa = 1./kappa_r;
const double kD = 4*Pi*sigma*D;
const double kD_r = 1.0/kD;
const double N_disso = 1 + ka*kD_r;



int status;
int iter = 0;
int max_iter = 100;
const gsl_root_fsolver_type *T;
gsl_root_fsolver *s;
double r = 0;
double r_expected = 6.26;
gsl_function F;
struct hx_params params = {  kappa, kappa1_r, kd, N_disso, gamma_b, gamma_u};
F.function = &hx;
F.params = &params;
T = gsl_root_fsolver_brent;
s = gsl_root_fsolver_alloc(T);
gsl_root_fsolver_set (s, &F, x_low, x_high);


       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           x_low = gsl_root_fsolver_x_lower (s);
           x_high = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (x_low, x_high, 0, 0.001);

         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fsolver_free (s);

    return r;

}




double find_roots_hy(phys_parameters_bimol2 physparams, double q_low, double q_high){


const double ka = physparams.ka;
const double D = physparams.D;
const double sigma = physparams.sigma;
const double b = physparams.b;
const double kd = physparams.kd;
const double gamma_b = physparams.gamma_b;
const double gamma_u = physparams.gamma_u;
const double sqrt_D = sqrt(D);
const double sqrt_D_r = 1.0/sqrt_D;
const double kappa = sqrt_D/(b-sigma);
const double kappa1_r = sigma*sqrt_D_r; 
const double kD = 4*Pi*sigma*D;
const double kD_r = 1.0/kD;
const double N_disso = 1 + ka*kD_r;



int status;
int iter = 0;
int max_iter = 100;
const gsl_root_fsolver_type *T;
gsl_root_fsolver *s;
double r = 0.0;


gsl_function F;
struct hy_params params = { kappa, kappa1_r, kd, N_disso, gamma_b, gamma_u};
F.function = &hy;
F.params = &params;
T = gsl_root_fsolver_brent;
s = gsl_root_fsolver_alloc(T);
gsl_root_fsolver_set (s, &F, q_low, q_high);


       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           q_low = gsl_root_fsolver_x_lower (s);
           q_high = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (q_low, q_high, 0, 0.0001);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fsolver_free (s);

       if (r < 10e-10){
        std::cout << __FILE__ << ": " <<__LINE__ <<  ": The first root almost zero!" << std::endl;
       }

    return r;

}

double find_roots_htan(phys_parameters_bimol2 physparams, double q_low, double q_high){


const double ka = physparams.ka;
const double D = physparams.D;
const double sigma = physparams.sigma;
const double b = physparams.b;
const double kd = physparams.kd;
const double gamma_b = physparams.gamma_b;
const double gamma_u = physparams.gamma_u;
const double sqrt_D = sqrt(D);
const double sqrt_D_r = 1.0/sqrt_D;
const double kappa = sqrt_D/(b-sigma);
const double kappa1_r = sigma*sqrt_D_r; 
const double kD = 4*Pi*sigma*D;
const double kD_r = 1.0/kD;
const double N_disso = 1 + ka*kD_r;



int status;
int iter = 0;
int max_iter = 100;
const gsl_root_fsolver_type *T;
gsl_root_fsolver *s;
double r = 0.0;


gsl_function F;
struct hy_params params = { kappa, kappa1_r, kd, N_disso, gamma_b, gamma_u};
F.function = &htan;
F.params = &params;
T = gsl_root_fsolver_brent;
s = gsl_root_fsolver_alloc(T);
gsl_root_fsolver_set (s, &F, q_low, q_high);


       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           q_low = gsl_root_fsolver_x_lower (s);
           q_high = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (q_low, q_high, 0, 0.001);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fsolver_free (s);

    return r;

}




double hx(double l, void *params){

  struct hx_params *p = (struct hx_params *) params;

  double result = 0.0;
  double x = 0.0;

  const double kappa1_r = p->kappa1_r;
  const double kappa = p->kappa;
  const double kd = p->kd;
  const double N_disso = p->N_disso;
  const double gamma_b = p->gamma_b;
  const double gamma_u = p->gamma_u;
  double kappa_sqrd = kappa*kappa;
  double kappa_cube = kappa_sqrd*kappa;

  x = 2*Pi*l;

  result = (kappa_cube*kappa1_r*x*x*x  + (kd -(gamma_u - gamma_b))*kappa*kappa1_r*x )*cosh(x) + (N_disso*(kappa_sqrd*x*x - (gamma_u-gamma_b))+ kd)*sinh(x);

 if ( fabs(result) > epsilon){
          return result;
     }
     else{
          return 0.0;
     }
}



double hy(double q, void *params){

struct hy_params *p  = (struct hy_params *) params;

double result = 0.0;
double y = 0.0;

const double kappa = p->kappa;
const double kappa1_r = p->kappa1_r;
const double kd = p->kd;
const double N_disso = p->N_disso;
const double gamma_b = p->gamma_b;
const double gamma_u = p->gamma_u;
double kappa_sqrd = kappa*kappa;
double kappa_cube = kappa_sqrd*kappa;

y = 2*Pi*q;

result = kappa1_r*(kappa_cube*y*y*y + (gamma_u-gamma_b -kd)*kappa*y)*cos(y) + ( N_disso*(kappa_sqrd*y*y + gamma_u-gamma_b)- kd)*sin(y);
 if ( fabs(result) > epsilon){
          return result;
     }
     else{
          return 0.0;
     }

}

double htan(double q, void *params){

struct hy_params *p  = (struct hy_params *) params;

double result = 0.0;
double y = 0.0;

const double kappa = p->kappa;
const double kappa1_r = p->kappa1_r;
const double kd = p->kd;
const double N_disso = p->N_disso;
const double gamma_b = p->gamma_b;
const double gamma_u = p->gamma_u;
double kappa_sqrd = kappa*kappa;
double kappa_cube = kappa_sqrd*kappa;

y = 2*Pi*q;

result = sin(y) + kappa*kappa1_r*y*cos(y) ;
 if ( fabs(result) > epsilon){
          return result;
     }
     else{
          return 0.0;
     }
}

  



double sketch_hx(double l, double kappa, double kappa1_r, double kd, double N_disso, double gamma_b, double gamma_u){

  double result = 0.0;
  double x = 0.0;
  
  double kappa_sqrd = kappa*kappa;
  double kappa_cube = kappa_sqrd*kappa;

  x = 2*Pi*l;
  
  result = (kappa_cube*kappa1_r*x*x*x  + (kd -(gamma_u - gamma_b))*kappa*kappa1_r*x )*cosh(x) + (N_disso*(kappa_sqrd*x*x - (gamma_u-gamma_b))+ kd)*sinh(x);
  
  return result;

}



double sketch_hy(double q, double kappa, double kappa1_r, double kd, double N_disso, double gamma_b, double gamma_u){

  double result = 0.0;
  double y = 0.0;
  
  double kappa_sqrd = kappa*kappa;
  double kappa_cube = kappa_sqrd*kappa;

  y = 2*Pi*q;
  
  result = kappa1_r*(kappa_cube*y*y*y + (gamma_u - gamma_b - kd)*kappa*y)* cos(y) + (N_disso*(kappa_sqrd*y*y + gamma_u - gamma_b) - kd)*sin(y);
  
  return result;

}













