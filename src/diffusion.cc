/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- Martin Schlather, Reinhard Furrer, Martin Kroll

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
 
#include "Basic_utils.h"  // must be before anything else
#ifdef DO_PARALLEL
#include <omp.h>
#endif
#include <R.h>
#include <R_ext/Lapack.h>
#include "RandomFieldsUtils.h"
#include "own.h"
#include "init_RandomFieldsUtils.h"
#include "kleinkram.h"

SEXP Udiffusion(SEXP SUSc, SEXP SUCo, SEXP Snevertried,
		SEXP Sa, SEXP Sabar, SEXP StWeight, SEXP Sq, SEXP Sdt,
		SEXP Srho, SEXP SrandSc, SEXP SrandCo, SEXP Sit, SEXP Sdummy,
		SEXP Sthreshold) {
  #define r_per_step 2
  int
    *nevertried = INTEGER(Snevertried),
    repN = length(SUSc),
    N = nrows(Snevertried),
    rep = ncols(Snevertried),
    it = INTEGER(Sit)[0]
   ;
  double
    *USc = REAL(SUSc),
    *UCo = REAL(SUCo),
    *a = REAL(Sa),
    *abar = REAL(Sabar),
    *tWeight = REAL(StWeight),
    *q = REAL(Sq),
    dt = REAL(Sdt)[0],
    *randSc = REAL(SrandSc) + it * repN,
    *randCo = REAL(SrandCo) + it * repN,
    *dummy = REAL(Sdummy),
    rho = REAL(Srho)[0],
    threshold = REAL(Sthreshold)[0]
    ;

  SEXP Ans;
  PROTECT(Ans = allocVector(INTSXP, rep));
  int *dN = INTEGER(Ans);
  GetRNGstate();
#ifdef DO_PARALLEL
#pragma omp parallel for 
#endif
  for (int r=0; r<rep; r++) {
    //
    //    printf("r=%d\n", r);
    double *usc = USc + N * r,
      *uco = UCo + N * r,
      *uscnew = dummy + N * r, // zwingend so gross wegen o m p !
      *RSc = randSc + N * r,
      *RCo = randCo + N * r,
      qr = q[r];
    int *never = nevertried + N * r,
       deltaN = 0;
    //    printf("ok %d\n", N);
    xA(usc, tWeight, N, N, uscnew);
    //    printf("ok\n");
    for (int i=0; i<N; i++) {
      // printf("i=%d\n", i);
      usc[i] += (rho * uscnew[i] + qr) * dt + RSc[i];
      uco[i] += RCo[i];
      if (never[i] && a[i] * usc[i] + abar[i] * uco[i] >= threshold) {
	never[i] = false;
	deltaN++;
      }
    }
    dN[r] = deltaN;
  }

  PutRNGstate();
  UNPROTECT(1);
  return(Ans);
}
