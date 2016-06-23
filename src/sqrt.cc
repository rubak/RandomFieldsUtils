/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- Martin Schlather, Reinhard Furrer

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
#include "Solve.h"
#include "own.h"


/*

## extrem wichter check -- funktioniert bei spam nicht
library(RandomFields, lib="~/TMP")
RFoptions(printlevel = 3, pch="", seed=999, use_spam = TRUE)
z = RFsimulate(RMspheric(), x, max_variab=10000, n=10000, spC=FALSE)
C = cov(t(z))
c = RFcovmatrix(RMspheric(), x)
print(summary(as.double(c - C))) ##//
stopifnot(max(abs(c-C)) < 0.05)

 */


int sqrtPosDef(double *M, int size,    // in out
	  solve_storage *pt     // in out
	       ){
  int err, sizeSq = size * size;
  //  InversionMethod Methods[SOLVE_METHODS] = { GLOBAL.solve.Methods[0], 
  //					     GLOBAL.solve.Methods[1] };
  //  GLOBAL.solve.Methods[0] = GLOBAL.solve.Methods[1] = 
 
  //int size2 = MIN(5, size);   printf("Mcall (%d x %d)=\n", size, size);for (int ii=0; ii<size2; ii++) {for (int jj=0; jj<size2; jj++) {printf("%e ", M[ii + jj * size]);}printf("\n");}

  if (GLOBAL.solve.sparse == True) 
    warning("package 'spam' is currently not used for simulation");
  usr_bool sparse = GLOBAL.solve.sparse;
  GLOBAL.solve.sparse = False;
  assert(pt != NULL);
  CMALLOC(result, sizeSq, double);
  err = doPosDef(M, size, true, NULL, 0, result, NULL, true, pt,
		 &(GLOBAL.solve));

  // printf("err = %d\n", err);BUG;

  GLOBAL.solve.sparse = sparse;
  return err;
}

static solve_param chol_param = chol_param_default;
SEXP CholPosDef(SEXP M) {
  return doPosDef(M, R_NilValue, R_NilValue, true, &chol_param);
}

int sqrtRHS(solve_storage *pt, double* RHS, double *result){

  // printf("%ld\n", (long) pt);  printf("SIZE %d\n", pt->size);
 
  assert(pt != NULL);
  int k = 0,
    size = pt->size;
  switch (pt->method) { 
  case direct_formula : 
  case Cholesky : {
    assert(pt->U != NULL);
    double *U = pt->result;
    //   for (int i=0; i<size; i++) {
    //  for(int j=0; j<size; j++) printf("%6.3f ", U[i+j*size]); printf("\n");} 

    for (int i=0; i<size; i++, k+=size) {
      double *Uk = U + k,
	dummy = 0.0;
      for (int j=0; j<=i; j++) dummy += RHS[j] * Uk[j];
      result[i] = (double) dummy; 
      //     result[i] = i; printf("nonsense %d\n", i);
    }
  }
    break;

  case SVD : {  
    assert(pt->U != NULL);
    double *U = pt->result;
    //  for (int i=0; i<size; i++) { for(int j=0; j<size; j++) printf("%6.3f ", U[i + j * size]); printf("\n");  }
  
    for (int i=0; i<size; i++){
      double dummy = 0.0;
      k = i;
      for (int j=0; j<size; j++, k+=size) dummy += U[k] * RHS[j];
      result[i] = (double) dummy; 
    }
  }
    break;

  case Sparse : {
    BUG; // SEE ALSO solve, sqrtOnly, tmp_delete !!
    int one = 1;
    F77_CALL(amuxmat)(&size, &size, &one, RHS, pt->DD, pt->lnz, 
		      pt->xja, pt->xlnz);
    //for (int i=0; i<size; i++) printf("%d %f %f\n", pt->invp[i], RHS[i], pt->DD[i]);
    for (int i=0; i<size; i++) result[i] = pt->DD[pt->invp[i]];

    // for (int i=0; i<size; i++) printf("%d %f %f\n",  pt->invp[i], RHS[i], result[i]);
  }
    break;

  case Diagonal : {  
    int  i, j,
      sizeP1 = size + 1;
    double *D = pt->result;
    for (i=j=0; j<size; j++, i+=sizeP1) result[j] = RHS[j] * D[i];
  }
    break;
  default : 
    //        printf("emthod %d %s %d\n", pt->method, InversionNames[pt->method], size);
    BUG;
  }
  
  return NOERROR;
}
