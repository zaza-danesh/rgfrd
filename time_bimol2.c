//
//  time_bimol2.c
//  
//
//  Created by Zahedeh Bashardanesh on 05/26/14.
//
//

//Z: the probability of quiting the interactive zone; -dS(t)/dt, this method is returning survival_rev retursn : S(t)
 
#include "time_bimol2.h"

 

double time_bimol2(double t, It_coeff_struct It_coeff, Qt_coeff_struct Qt_coeff){

   double result = 0.0;
   double It = 0.0;
   double Qt = 0.0;

   It = It_bimol2(t, It_coeff);
   Qt = Qt_bimol2(t, Qt_coeff);
   result = 1.0 - (It + Qt);
   
   return result;
}




