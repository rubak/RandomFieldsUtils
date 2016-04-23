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

const char * solve[solveN] = SOLVE_NAMES;	

int sqrt_(double *M, int size,    // in out
	  solve_storage *pt,      // in out
	  solve_param *Sp, int PL // in
	  ){
  /*

    M: positiv definite square matrix (symmetry is not checked) of size x size;
       NOTE THAT THE CONTENTS OF M IS DESTROYED, the sqrt of M is returned here
    pt  : working space. If NULL, internal working spaces are used.
          A non-NULL value gives an advantage only if sqrt_ is called
          repeated. Then
            solve_storage *pt = (solve_storage*) malloc(sizeof(solve_storage);
            solve_NULL(pt);
          prepares pt. Deletion is done only at the very end by
            solve_DELETE(pt);
          In meanwhile the working space is managed by sqrt__;
     Sp : parameters. If NULL, the behaviour is as described in
          the R manual for sqrt_.
          The parameters are described in the package 'RandomFields'
          under ?RFoptions
     PL : printing level; 1 is standard
	  
  */

  assert(SOLVE_METHODS >= 2);
  solve_param
    sp_standard = solve_param_default,
    *sp = Sp == NULL ? &sp_standard : Sp;
  int
    err = NOERROR,
    spam_zaehler = 0,
    nnzA = 0,
    sizeSq = size * size,
    sizeP1 = size + 1;
  usr_bool sparse = sp->sparse;
  double  
    svdtol = sp->svd_tol,
    spam_tol = sp->spam_tol;
  bool diag  = false;
 
  pt->method = NoFurtherInversionMethod;
  pt->size = size;

  if (sparse == Nan && (sparse = (usr_bool) (size > sp->spam_min_n))) {
    double mean_diag = 0.0;
    for (int i=0; i<sizeSq; i += sizeP1) mean_diag += M[i];
    mean_diag /= (double) size;
    spam_tol *= mean_diag;

    bool random_sample = sizeSq >= sp->spam_sample_n * 3;
    if (random_sample) {
      double 
	thr = sp->spam_sample_n * (1.0 - sp->spam_min_p);
      int	
	threshold = (int) (thr + sqrt(thr) * 3),
	notZero = 0;
      for (int i=0; i<sp->spam_sample_n; i++) {
	if ((notZero += fabs(M[(i * sp->spam_factor) % sizeSq]) > 
	     spam_tol) >= threshold){
	  sparse = False;
	  break;
	}
      }
      if (PL >= PL_FCTN_DETAILS) 
	PRINTF("random sampling: sparse=%d\n", 
	       sparse == Nan ? NA_INTEGER : (int) sparse);
    }
    if (!random_sample || sparse == True) {
      int diag_nnzA = 0;
      for (int i=0; i<size; i++) {
	int end = i * sizeP1;
	long j;
	for (j=i * size; j<end; nnzA += fabs(M[j++]) >= spam_tol);
	diag_nnzA += fabs(M[j++]) > spam_tol;
	end = (i+1) * size;
      }
      diag = (nnzA == 0);
      nnzA *= 2;
      nnzA += diag_nnzA;
      sparse = (usr_bool) (nnzA <= sizeSq * (1.0 - sp->spam_min_p));
      spam_zaehler = nnzA + 1;
      if (PL >= PL_DETAILSUSER) {
	if (diag) PRINTF("diagonal matrix detected\n");
 	else if (sparse == True) 
	  PRINTF("sparse matrix detected (%3.2f%% zeros)\n", 
		 100.0 * (1.0 - nnzA / (double) sizeSq));
	else PRINTF("full matrix detected (%3.2f%% nonzeros)\n", 
		    100.0 * nnzA / (double) sizeSq);
      }
    }
  } else {
    diag = true;
    for (int i=0; i<size && diag; i++) {
      int end = i * sizeP1;
      long j;
      for (j=i * size; j<end; j++)
	if (fabs(M[j]) > spam_tol) {
	  diag = false;
	  break;
	}
      j++;
      end = (i+1) * size;
    }
  }
  
  if (pt == NULL) BUG;

  if (diag) {
    pt->method = Diagonal;
    CMALLOC(D, size, double); 
    for (int i=0, j=0; i<sizeSq; i += sizeP1, j++)
      D[j] = M[i] > 0.0 ? sqrt(M[i]) : 0.0;
    err = NOERROR;
    goto ErrorHandling;
  }
  
  InversionMethod *Meth;
  if (sparse == True || sp->Methods[0] >= NoFurtherInversionMethod) {
    Meth = pt->newMethods;
    if (sparse == True) {
      pt->newMethods[0] = Sparse;
      pt->newMethods[1] = 
	sp->Methods[0] < NoFurtherInversionMethod && sp->Methods[0] != Sparse
	? sp->Methods[0] 
	: Cholesky;  
      if (SOLVE_METHODS > 2) {
	pt->newMethods[2] = 
	  sp->Methods[0]<NoFurtherInversionMethod && sp->Methods[0] != Sparse &&
	  sp->Methods[1]<NoFurtherInversionMethod && sp->Methods[1] != Sparse
	  ? sp->Methods[1] 
	  : SVD;
      }
      // pt->newMethods[1] = Sparse;
    } else {
      pt->newMethods[0] = Cholesky;  
      pt->newMethods[1] =  SVD;
      if (SOLVE_METHODS > 2) pt->newMethods[2] = SVD;
    }
  } else Meth = sp->Methods;

  for (int m=0; m<SOLVE_METHODS && (m==0 || Meth[m] != Meth[m-1]); m++) {
    pt->method = Meth[m];
    // 
    if (pt->method == NoFurtherInversionMethod ||
	pt->method == NoInversionMethod) break;
   if (PL>=PL_STRUCTURE) { 
      PRINTF("method to calculate square root : %s\n", 
	     InversionNames[pt->method]);
    }
 
    // print("m=%d Method %d chol=%d QR=%d SVD=%d sparse=%d\n", m, method, Cholesky, QR, SVD, Sparse);
   
    
    switch(pt->method) {
 
    case Cholesky : { // cholesky       
       if (size > sp->max_chol)
	 CERR("matrix too large for Cholesky decomposition.");
       CMALLOC(U, sizeSq, double);
       memcpy(U, M, sizeSq * sizeof(double));
      
       F77_CALL(dpotrf)("U", &size, U, &size, &err);  
	 
       /*
	 if (false)  {
	 double *sq  = (double *) MALLOC(sizeof(double) * size * size);
	 AtA(U, size, size, sq);
	 
	 long i,j;
	 PRINTF("AtA \n");
	   for (int i=0; i<locpts * vdim; i++) {
	   for (int j=0; j<locpts * vdim; j++) {
	   PRINTF("%+2.2f ", Cov[i  + locpts * vdim * j]);
	   }
	   PRINTF("\n");
	   }
	   //  assert(false);
	   }
       */
       
	 if (err == NOERROR) {
	   F77_CALL(dpotri)("U", &size, M, &size, &err);  
	   long  
	     i2=0;
	   for (int i=0; i<size; i++, i2+=size + 1) {	
	     for (long i3 = i2 + 1, j = i2 + size; j<sizeSq; j+=size)
	       M[i3++] = M[j];
	   }
	   
	 }
	  
	 if (err != NOERROR) {	
	   CERR1("Cholesky decomposition failed with err=%d of 'dpotr*' in 'sqrt_'. Probably matrix not strictly positive definite.\n", err);
	 }
	 
       if (PL >= PL_DETAILSUSER) PRINTF("Cholesky decomposition successful\n");
    }
      break;

    case SVD : {// SVD : M = U D VT

      if (size > sp->max_svd) CERR("matrix too large for SVD decomposition.");
      double  optim_lwork;
      int k = 0,  
	size8 = size * 8,
	lwork = -1;
  
      CMALLOC(DD, sizeSq, double);
      memcpy(DD, M, sizeSq * sizeof(double));
      
      CMALLOC(VT, sizeSq, double);
      CMALLOC(U, sizeSq, double);
      CMALLOC(D, size, double); 
      CMALLOC(iwork, size8, int);
  
      F77_CALL(dgesdd)("A", &size, &size, DD, &size, D, U, &size, VT, &size, 
		       &optim_lwork, &lwork, iwork, &err);
      if (err != NOERROR) {
	CERR1("'dgesdd' failed with err=%d\n", err);
	break;
      }

      lwork = (int) optim_lwork;
      CMALLOC(work, lwork, double);

      F77_CALL(dgesdd)("A",  &size,  &size, DD, &size, D, U, &size, VT, &size, 
		       work, &lwork, iwork, &err);
      
      if (err != NOERROR || ISNAN(D[0])) {
	if (PL>PL_ERRORS) PRINTF("Error code F77_CALL(dgesdd) = %d\n", err); 
	CERR1("'dgesdd' failed with err=%d\n", err);
	break;
      }
         
      /* calculate SQRT of covariance matrix */
      for (int j=0; j<size; j++) {
	double dummy;
        dummy = fabs(D[j]) < svdtol ? 0.0 : sqrt(D[j]);
	for (int i=0; i<size; i++) {
	  U[k++] *= dummy;
	}
      }

      /* check SVD */
      if (svdtol >= 0.0) {
	for (int i=0; i<size; i++) {
	  for (k=i; k<size; k++) {
	    double sum = 0.0;
	    for (int j=0; j<sizeSq; j+=size) sum += U[i+j] * U[k+j];
	    
	    if (fabs(M[i * size + k] - sum) > svdtol) {
	      if (PL > PL_ERRORS) {
		PRINTF("difference %e at (%d,%d) between the value (%e) of the covariance matrix and the square of its root (%e).\n", 
		       M[i * size +k] - sum, i, k, M[i* size +k ], 
		       sum);
	      }
	      FERR3("required precision not attained  (%e !> %e): probably invalid model. See also '%s'\n", fabs(M[i * size + k] - sum), svdtol,
		    solve[SOLVE_SVD_TOL]);
	      err=ERRORM;
	      break;
	    } //else printf("ok (%d,%d) %f %f\n", i, k, M[i* size +k ],sum);
	  }
	  if (err != NOERROR) break;		
	}
	if (err != NOERROR) break;		
     }
     
      if (PL >=  PL_DETAILSUSER) PRINTF("svd successful\n");
      break;
    }

    case Sparse : {// sparse matrix

     int 
       doperm = sp->pivot,
       halfsq = size * (size + 1) / 2,
       nnzcolindices = 0,
       nnzR = 0,
       cache = 512, // to do: CPU cache size
       nnzcfact[3] = { 5, 1, 5 }, 
       nnzRfact[3] = { 5, 1, 2 };
      double
	cholincrease_nnzcol = 1.25,
        cholincrease_nnzR = 1.25;

      CMALLOC(pivot, size, int);
      if (!doperm) for (int i=0; i<size; i++) pivot[i] = i + 1;

      if (spam_zaehler == 0) {
	for (int i=0; i<sizeSq; i++) nnzA += fabs(M[i]) >= spam_tol;
	spam_zaehler = nnzA + 1; // falls nur aus Nullen bestehend
      }
      
      CMALLOC(xlnz, sizeP1, int);
      CMALLOC(snode, size, int);
      CMALLOC(xsuper, sizeP1, int);
      CMALLOC(xlindx, sizeP1, int);
      CMALLOC(invp, size, int);
      CMALLOC(cols, spam_zaehler, int);
      CMALLOC(rows, sizeP1, int);
   
      int nDD = spam_zaehler;
      if (nDD < size) nDD = size;
      CMALLOC(DD, nDD, double);
     // prepare spam

      F77_CALL(spamdnscsr)(&size, &size, M, &size, DD, cols, rows, 
			   &spam_tol); // create spam object
     
     // calculate spam_cholesky
      err = 4; // to get into the while loop
      while (err == 4 || err == 5) {
	if (nnzcolindices == 0) {
	  double rel = nnzA / (double) size;
	  if (rel < 5) {
	    nnzcolindices = (int) ceil(nnzA * (1.05 * rel - 3.8));
	    if (nnzcolindices < 1000) nnzcolindices = 1000;
	  } else {
	    nnzcolindices = nnzA;
	  }
	  nnzcolindices *= nnzcfact[doperm];
	  if (nnzcolindices < nnzA) nnzcolindices = nnzA;
	} else if (err == 5) {
	  int tmp = (int) ceil(nnzcolindices * cholincrease_nnzcol);
	  if (PL > PL_RECURSIVE) 
	    PRINTF("Increased 'nnzcolindices' with 'NgPeyton' method\n(currently set to %d from %d)", tmp, nnzR);
	  nnzcolindices = tmp;
	}
	if (nnzcolindices < pt->lindx_n) nnzcolindices = pt->lindx_n;
	
	if (nnzR == 0) {
	  double u = floor(.4 * pow(nnzA, 1.2));
	  u = u < 4 * nnzA ? 4 * nnzA : ceil(u);
	  nnzR = (int) u * nnzRfact[doperm];
	} else if (err == 4) {
	  int tmp = (int) ceil(nnzR * cholincrease_nnzR);
	  if (PL > PL_RECURSIVE) 
	    PRINTF("Increased 'nnzR' with 'NgPeyton' method\n(currently set to %d from %d)", tmp, nnzR);
	  nnzR = tmp;
	}
	if (nnzR < pt->lnz_n) nnzR = pt->lnz_n;
	else if (nnzR > halfsq) nnzR = halfsq;	
	
	CMALLOC(lindx, nnzcolindices, int);	
	CMALLOC(lnz, nnzR, double);
	
	F77_CALL(cholstepwise)(&size, &nnzA, DD, cols, rows, &doperm,
			       invp, pivot, 
			       &(pt->nnzlindx), &nnzcolindices, 
			       lindx, // 
			       xlindx,// 
			       &(pt->nsuper), // length of lindx
			       &nnzR,  // physical length of lindx
			       lnz,   // output:result
			       xlnz,  // cols of lnz
			       snode,  // supernode membership ??
			       xsuper, // supernode partioning
			       &cache, // cache size of the CPU
			       &err
			       );       
	
	if (err != NOERROR) {
	  CERR1("'cholstepwise' failed with err=%d\n", err);
	  break;
	}	 
      } // while
      
      if (err != NOERROR) CERR("'spam' failed");     
      if (PL >=  PL_DETAILSUSER) PRINTF("'spam' successful\n");

      nnzR = xlnz[size] - 1;
      CMALLOC(xja, nnzR, int);
      F77_CALL(calcja)(&size, &(pt->nsuper), pt->xsuper, 
		       pt->lindx, pt->xlindx, pt->xlnz, 
		       xja);
      for (int i=0; i<size; invp[i++]--);
      
      break;
    }
   
    default :
      GERR("unknown method for 'sqrt_'");
      
    } // switch
    if (err==NOERROR) break;
  } // for m

 ErrorHandling:
   return err; //  -method;
}
  

int sqrt_RHS_(solve_storage *pt, double* RHS, double *res){
  int k = 0,
    size = pt->size;
  switch (pt->method) {
  case Cholesky : {
    double *U = pt->U;
    //   for (int i=0; i<size; i++) {
    //    for(int j=0; j<size; j++) printf("%6.3f ", U[i+j*size]); printf("\n");}

    for (int i=0; i<size; i++, k+=size) {
      double *Uk = U + k,
	dummy = 0.0;
      for (int j=0; j<=i; j++) dummy += RHS[j] * Uk[j];
      res[i] = (double) dummy; 
    }
  }
    break;

  case SVD : {  
    double *U = pt->U;
    //  for (int i=0; i<size; i++) { for(int j=0; j<size; j++) printf("%6.3f ", U[i + j * size]); printf("\n");  }
  
    for (int i=0; i<size; i++){
      double dummy = 0.0;
      k = i;
      for (int j=0; j<size; j++, k+=size) dummy += U[k] * RHS[j];
      res[i] = (double) dummy; 
    }
  }
    break;

  case Sparse : {
    int one = 1;
    F77_CALL(amuxmat)(&size, &size, &one, RHS, pt->DD, pt->lnz, 
		      pt->xja, pt->xlnz);
    //for (int i=0; i<size; i++) printf("%d %f %f\n", pt->invp[i], RHS[i], pt->DD[i]);
    for (int i=0; i<size; i++) res[i] = pt->DD[pt->invp[i]];

    // for (int i=0; i<size; i++) printf("%d %f %f\n",  pt->invp[i], RHS[i], res[i]);
  }
    break;

  case Diagonal : {
    double *D = pt->D;
    for (int i=0; i<size; i++) res[i] = RHS[i] * D[i];
  }
    break;
  default : 
    // printf("emthod %d\n", pt->method);
    BUG;
  }
  return NOERROR;
}
