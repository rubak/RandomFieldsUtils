/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/
#include <R_ext/Lapack.h>
#include "RandomFieldsUtils.h"


double struve_(double x, double nu, double factor_Sign, bool expscaled)
{ 
 if ((x == 0.0) && (nu>-1.0)) return 0.0;
  if (x <= 0.0) return RF_NA; // not programmed yet
  double exp_dummy,
     dummy = 0.0, 
     logx = 2.0 * log(0.5 * x), 
     x1 = 1.5, 
     x2 = nu + 1.5,
     value = 1.0, 
     fsign = factor_Sign,
    epsilon=1e-20;

   do {
     dummy += logx - log(x1) - log(fabs(x2));
     exp_dummy = exp(dummy);
     value += (1 - 2 * (x2 < 0))  * fsign * exp_dummy;
     //  printf("%f %f %f %f\n", value, fsign, x1, x2);
     x1 += 1.0;
     x2 += 1.0;
     fsign = factor_Sign * fsign; 
   } while (exp_dummy > fabs(value) * epsilon);

 

   x1 = 1.5;
   x2 = nu + 1.5;
   if (x2 > 0.0) { 
     dummy = (nu + 1.0) * 0.5 * logx - lgammafn(x1) - lgammafn(x2);
     if (expscaled) dummy -= x;
     value *= exp(dummy);
   } else {
     //if ( (double) ((int) (x1-0.5)) != x1-0.5 ) return RF_NA;
     value *= pow(0.5 * x, nu + 1.0) / (gammafn(x1) * gammafn(x2));
     if (expscaled) value *= exp(-x);
   }

  return value;
}

//void StruveH(double *x, double *nu) {*x=struve(*x, *nu, -1.0, false);}
//void StruveL(double *x, double *nu, int * expScaled) {
//  *x=struve(*x, *nu, 1.0, (bool) *expScaled);
//}
double StruveH(double x, double nu) {return struve_(x, nu, -1.0, false);}
double StruveL(double x, double nu, bool expScaled) {
 return struve_(x, nu, 1.0, expScaled);
}


SEXP struve(SEXP X, SEXP Nu, SEXP Factor_Sign, SEXP Expscaled) {
  int i,    
    lenx = length(X),
    lennu = length(Nu),
    len = lenx;  
  if (len < lennu) len = lennu;
  SEXP Result;
  PROTECT(Result = allocVector(REALSXP, len));
  double *x = REAL(X),
    *nu  = REAL(Nu),
    factor_sign = REAL(Factor_Sign)[0],
    *result = REAL(Result);
  bool expscaled = LOGICAL(Expscaled)[0];
  for (i=0; i<len; i++)
    result[i] = struve_(x[i % lenx], nu[i % lennu], factor_sign, expscaled);
 
  UNPROTECT(1);
  return Result;
}



//void I0ML0(double *x, int *n) {
//  int i;
//  for (i=0; i<*n; i++) x[i] = I0mL0(x[i]);
//} 


double I0mL0(double x){
  /* Bessel I_0 - Struve L_0 for non-negative arguments x */
  /* see J. MacLeod, Chebyshev expansions for modified {S}truve and 
     related functions, Mathematics of Computation, 60, 735-747, 1993 */
  static double g2[24] = {
	0.52468736791485599138e0,
	-0.35612460699650586196e0,
	0.20487202864009927687e0,
	-0.10418640520402693629e0,
	0.4634211095548429228e-1,
	-0.1790587192403498630e-1,
	0.597968695481143177e-2,
	-0.171777547693565429e-2,
	0.42204654469171422e-3,
	-0.8796178522094125e-4,
	0.1535434234869223e-4,
	-0.219780769584743e-5,
	0.24820683936666e-6,
	-0.2032706035607e-7,
	0.90984198421e-9,
	0.2561793929e-10,
	-0.710609790e-11,
	0.32716960e-12,
	0.2300215e-13,
	-0.292109e-14,
	-0.3566e-16,
	0.1832e-16,
	-0.10e-18,
	-0.11e-18
  };
  static double g3[24] = {
	2.00326510241160643125e0,
	0.195206851576492081e-2,
	0.38239523569908328e-3,
	0.7534280817054436e-4,
	0.1495957655897078e-4,
	0.299940531210557e-5,
	0.60769604822459e-6,
	0.12399495544506e-6,
	0.2523262552649e-7,
	0.504634857332e-8,
	0.97913236230e-9,
	0.18389115241e-9,
	0.3376309278e-10,
	0.611179703e-11,
	0.108472972e-11,
	0.18861271e-12,
	0.3280345e-13,
	0.565647e-14,
	0.93300e-15,
	0.15881e-15,
	0.2791e-16,
	0.389e-17,
	0.70e-18,
	0.16e-18
  };
  double r, x2, ac;
  int i;
  
  if (x < 0.0) {return RF_NA;}
  if (x < 16.0) {
    r = 0.5 * g2[0];
    ac = acos((6.0 * x - 40.0) / (x + 40.0));
    for (i=1; i<24; i++) {
      r += g2[i] * cos(i * ac);
    }
  } else {
    r = 0.5 * g3[0];
    x2 = x * x;
    ac = acos((800.0 - x2) / (288.0 + x2));
    for (i=1; i<24; i++) {
      r += g3[i] * cos(i * ac);
    }
    r *= T_PI /* 2/pi */ / x;
  }
  return r;
}


SEXP I0ML0(SEXP X) {
  SEXP Result;
  PROTECT(Result = allocVector(REALSXP, length(X)));
  double *x = REAL(X),
    *result = REAL(Result);
  int i,    
    lenx = length(X);  
  for (i=0; i<lenx; i++) result[i] = I0mL0(x[i]);

  UNPROTECT(1);
  return Result;
}
