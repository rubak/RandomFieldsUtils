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

void solve_DELETE0(solve_storage *x) {
    FREE(x->iwork);
    FREE(x->ipiv);

    FREE(x->pivot);
    FREE(x->xlnz);
    FREE(x->snode);
    FREE(x->xsuper);
    FREE(x->xlindx);
    FREE(x->invp);
    FREE(x->cols);
    FREE(x->rows);
    FREE(x->lindx);
    //    FREE(x->);
   
    // double *
    FREE(x->SICH);
    FREE(x->MM);
    FREE(x->workspaceD);
    FREE(x->workspaceU);

    FREE(x->VT);
    FREE(x->work);
    FREE(x->U);
    FREE(x->D);

    FREE(x->workLU);
  
    FREE(x->lnz); 
    FREE(x->DD);
    FREE(x->RHS);
    FREE(x->w3);
}

void solve_DELETE(solve_storage **S) {
  solve_storage *x = *S;
  if (x!=NULL) {
    solve_DELETE0(*S);
    UNCONDFREE(*S);
  }
}
void solve_NULL(solve_storage* x) {
  if (x == NULL) return;
    x->iwork_n = x->ipiv_n = 
     x->pivot_n = x->xlnz_n = x->snode_n = x->xsuper_n = x->xlindx_n = 
      x->invp_n = x->cols_n = x->rows_n =x->lindx_n = 
      //
    x->SICH_n = x->workspaceD_n = x->workspaceU_n = 
    x->VT_n = x->work_n = x->w2_n = x->U_n = x->D_n = x->workLU_n =
      x->lnz_n =  x->DD_n = x->RHS_n = x->w3_n =
    0;
  
    x->iwork = x->ipiv = 
    x->pivot = x->xlnz = x->snode = x->xsuper = x->xlindx = 
    x->invp = x->cols = x->rows = x->lindx = 
    NULL;

    x->SICH = x->MM = x->workspaceD = x->workspaceU = 
    x->VT = x->work = x->w2 = x->U = x->D = x->workLU =
    x->lnz = x->DD = x->RHS = x->w3 =
    NULL;
}


#define CMALLOC(WHICH, N, TYPE)				 \
  if (pt->WHICH##_n < (N)) {					 \
    FREE(pt->WHICH);						 \
    pt->WHICH##_n = N;					      	 \
      if ((pt->WHICH = (TYPE *) CALLOC(N, sizeof(TYPE))) == NULL)	\
      return ERRORMEMORYALLOCATION;				 \
  } else {							 \
    assert( ((N) > 0 && pt->WHICH != NULL) || (N) == 0);	 \
  }								 \
  TYPE VARIABLE_IS_NOT_USED *WHICH = pt->WHICH


int solvePosDef_(double *M, int size, bool posdef,
		double *rhs, int rhs_cols,
		double *logdet, 
		solve_storage *pt, solve_param *sp, int PL
		){

  // http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F01/f01intro.xml#
  assert(SOLVE_METHODS >= 2);

  int m, i,
    err = NOERROR,
    method = -1,
    spam_zaehler = 0,
    nnzA = 0,
    sizeSq = size * size,
    sizeP1 = size + 1;
  long j;
  double 
    sparse = sp->sparse,
    spam_tol = sp->spam_tol;

  if (ISNAN(sparse) || ISNA(sparse)) {
    if ((sparse = size > sp->spam_min_n)) {
 	bool random_sample = sizeSq >= sp->spam_sample_n * 3;
      if (random_sample) {
	double 
	  thr = sp->spam_sample_n * (1.0 - sp->spam_min_p);
	int	
	  threshold = (int) (thr + sqrt(thr) * 3),
	  notZero = 0;
	for (i=0; i<sp->spam_sample_n; i++) {
	  if ((notZero += fabs(M[(i * sp->spam_factor) % sizeSq]) > 
	       spam_tol) >= threshold){
	    sparse = false;
	    break;
	  }
	}
      }
      if (!random_sample || sparse) {
	for (i=0; i<sizeSq; i++) nnzA += fabs(M[i]) > spam_tol;
	sparse = nnzA <= sizeSq * (1.0 - sp->spam_min_p);
	spam_zaehler = nnzA + 1;
      }
    }
  }
  
 
  //printf("sparse = %f size=%d min=%d sp->sparse=%f\n", sparse, size, sp->spam_min_n, sp->sparse);

  int *Meth;
  if (sparse || sp->Methods[0] < 0) {
    Meth = pt->newMethods;
    if (sparse) {
      pt->newMethods[0] = pt->newMethods[1] = Sparse;
    } else {
      pt->newMethods[0] = posdef ? Cholesky : LU;  
      pt->newMethods[1] = LU;
      if (SOLVE_METHODS > 2) pt->newMethods[2] = LU;
    }
  } else Meth = sp->Methods;

  //  PL = 10;

  if (!posdef && Meth[0] != LU && Meth[0] != SVD) 
    return ERRORNOTPROGRAMMEDYET;
 
  if (Meth[1] != Meth[0] && rhs_cols==0) { // at least two different Methods in the list
    CMALLOC(SICH, sizeSq, double);
    memcpy(SICH, M, sizeSq);
  }
  double *SICH = pt->SICH;

   for (m=0; m<SOLVE_METHODS && (m==0 || Meth[m] != Meth[m-1]); m++) {
    method = Meth[m];
     
    //   for (int u=0; u<sizeSq; u++) printf("%f ", M[u]);
    // printf("m=%d %d\n", m, size);

    if (method<0) break;
    if (m > 0 && rhs_cols == 0) {
      memcpy(M, SICH, sizeSq);
    }

    switch(method) {
 
     case Cholesky : // cholesky
       if (!posdef) CERR("Cholesky needs positive derfinite matrix");
       if (rhs_cols == 0) {
	 F77_CALL(dpotrf)("U", &size, M, &size, &err);  
	 if (err == NOERROR) {
	   if (logdet != NULL) {
	     for (*logdet=0.0, i=0; i < sizeSq; i+=sizeP1) *logdet += log(M[i]);
	     *logdet *= 2;
	   }
	 }	 
	 F77_CALL(dpotri)("U", &size, M, &size, &err);  
	 if (rhs_cols == 0) {
	   long  i2, i3;
	   for (i2=i=0; i<size; i++, i2+=size + 1) {	
	     for (i3 = i2 + 1, j = i2 + size; j<sizeSq; j+=size) M[i3++] = M[j];
	   }
	 }
       } else {
	 CMALLOC(MM, sizeSq, double);
	 MEMCOPY(MM, M, sizeSq * sizeof(double));
 	 F77_CALL(dpotrf)("U", &size, MM, &size, &err);  
	 if (err == NOERROR) {
	   if (logdet != NULL) {
	     for (*logdet=0.0, i=0; i < sizeSq; i+=sizeP1)
	       *logdet += log(MM[i]);
	     *logdet *= 2;
	   }
	   F77_CALL(dpotrs)("U", &size, &rhs_cols, MM, &size, rhs, &size, &err);
	 }
       }
       
       if (err != NOERROR) {	
	 if (PL >= PL_COV_STRUCTURE) 
	   PRINTF("Cholesky decomposition failed; matrix does not seem to be strictly positive definite (err# = %d)\n", err);
	 CERR1("'dpotr*' failed with err=%d\n", err);
       }

       if (PL >=  PL_STRUCTURE) PRINTF("Cholesky decomposition successful\n");
      
      break;

    case QR : {// QR returns transposed of the inverse !!
      if (rhs_cols > 0 || logdet != NULL) {
	err = ERRORFAILED;
	continue;
      }

      err = ERRORNOTPROGRAMMEDYET; /// to do: clarify transposed !
      continue;

      CMALLOC(workspaceD, size, double);
      CMALLOC(workspaceU, size, double);

      // PseudoInverse:
      //    F77_CALL(f01blf)(&size, &size, &(GLOBAL.general.matrixtolerance),
       //		       M, &size, aijmax, &irank, inc, workspaceD, 
      //		       workspaceU, &size, workspaceDU, &err);
      F77_CALL(dgeqrf)(&size, &size,
		       M, &size, // aijmax, &irank, inc, workspaceD, 
		       workspaceU, workspaceD,  &size, &err);     
      
      if (err != NOERROR) {	
	CERR1("'dgeqrf' failed with err=%d\n", err);
      }
      if (PL >=  PL_STRUCTURE) PRINTF("Cholesky decomposition successful\n");
      break;
    }

    case SVD : {// SVD : M = U D VT
      double  optim_lwork;
      int k,  
	size8 = size * 8,
	lwork = -1;
      double *M0 = M;
      CMALLOC(MM, rhs_cols == 0 ? 0 : sizeSq, double);
      if (rhs_cols > 0) {
	MEMCOPY(MM, M, sizeSq * sizeof(double));
	M0 = MM;
      }

      CMALLOC(VT, sizeSq, double);
      CMALLOC(U, sizeSq, double);
      CMALLOC(D, size, double); 
      CMALLOC(iwork, size8, int);
  
      F77_CALL(dgesdd)("A", &size, &size, M0, &size, D, U, &size, VT, &size, 
		       &optim_lwork, &lwork, iwork, &err);
      if (err != NOERROR) {
	CERR1("'dgesdd' failed with err=%d\n", err);
	break;
      }

      lwork = (int) optim_lwork;
      CMALLOC(work, lwork, double);

      F77_CALL(dgesdd)("A",  &size,  &size, M0, &size, D, U, &size, VT, &size, 
		       work, &lwork, iwork, &err);
      
      if (err!=NOERROR || ISNAN(D[0])) {
	if (PL>PL_ERRORS) PRINTF("Error code F77_CALL(dgesdd) = %d\n", err); 
	CERR1("'dgesdd' failed with err=%d\n", err);
	break;
      }

      // calculate determinant 
      if (logdet != NULL) {
	for (*logdet=0.0, i = 0; i < size; *logdet += log(D[i++]));
      }


      //      for (k=0; k<sizeSq; k++) printf("%f ", U[k]); printf("\nVT=");
      //      for (k=0; k<sizeSq; k++) printf("%f ", VT[k]); printf("\n D=");
      //      for (k=0; k<size; k++) printf("%f ", D[k]); printf("\ninvD=");
	  

      double svd_tol = sp->svd_tol;
      for (j=0; j<size; j++) D[j] = fabs(D[j]) < svd_tol ? 0.0 : 1.0 / D[j];

      if (rhs_cols > 0) {
	int tot = size * rhs_cols;
	CMALLOC(w2, tot, double);	
	matmulttransposed(U, rhs, w2, size, size, rhs_cols);
	for (k=0; k<tot; )
	  for (i=0; i<size; i++) w2[k++] *= D[i];
	// for (k=0; k<rhs_cols * size; k++) printf("%f ", w2[k]); printf("\n");
	matmulttransposed(VT, w2, rhs, size, size, rhs_cols);
      } else {
	// calculate inverse of covariance matrix
	for (k=0, j=0; j<size; j++) {
	  double dummy = D[j];
	  for (i=0; i<size; i++) U[k++] *= dummy;
	}
	matmult_tt(U, VT, M, size, size, size); // V * U^T
      } 

      if (PL >=  PL_STRUCTURE) PRINTF("svd successful\n");
      break;
    }


    case LU : {// LU
      double *M0 = M;
      CMALLOC(MM, rhs_cols == 0 ? 0 : sizeSq, double);
      if (rhs_cols > 0) {
	MEMCOPY(MM, M, sizeSq * sizeof(double));
	M0 = MM;
      }
      
      CMALLOC(ipiv, size, int);		    
      F77_CALL(dgetrf)(&size, &size, M0, &size, ipiv, &err);
      if (err != NOERROR) {
	CERR1("'dgetrf' (LU) failed with err=%d\n", err);
      }
 
      if (logdet != NULL)
	for (*logdet=0.0, i = 0; i < sizeSq; i += sizeP1) *logdet += log(M0[i]);
      
      if (rhs_cols > 0) {
	F77_CALL(dgetrs)("N", &size, &rhs_cols, M0, &size, ipiv, 
			 rhs, &size, &err);
	if (err != NOERROR) {	
	  CERR1("'dgetrs' (LU) failed with err=%d\n", err);
	}
      } else {
	int lwork = -1;
	double dummy;
	F77_CALL(dgetri)(&size, M0, &size, ipiv, &dummy, &lwork, &err);	
	if (err != NOERROR) {	
	  CERR1("'dgetri' (LU) failed with err=%d\n", err);
	}
 	lwork = (int) dummy;
	CMALLOC(workLU, lwork, double);
	F77_CALL(dgetri)(&size, M0, &size, ipiv, workLU, &lwork, &err);	//M0=M
 	if (err !=  NOERROR) {	
	  CERR1("'dgetri' (LU) failed with err=%d.\n", err);
	}
      }
      
      if (PL >=  PL_STRUCTURE) PRINTF("LU decomposition successful\n");
      break;
    }
     
    case Sparse : {// sparse matrix
     int nnzlindx, nsuper, RHS_COLS, 	
	doperm = sp->pivot,
	halfsq = size * (size + 1) / 2,
	nnzcolindices = 0,
	nnzR = 0,
	cache = 512; // to do: CPU cache size
      double
	nnzcfact[3] = { 5.0, 1.0, 5.0 }, 
	nnzRfact[3] = { 5.0, 1.0, 2.0 },
        cholincrease_nnzcol = 1.25,
        cholincrease_nnzR = 1.25;

      CMALLOC(pivot, size, int);
      if (!doperm) for (i=0; i<size; i++) pivot[i] = i + 1;

      if (spam_zaehler == 0) {
	for (i=0; i<sizeSq; i++) nnzA += fabs(M[i]) > spam_tol;
	spam_zaehler = nnzA + 1; // falls nur aus Nullen bestehend
      }
      
      CMALLOC(xlnz, sizeP1, int);
      CMALLOC(snode, size, int);
      CMALLOC(xsuper, sizeP1, int);
      CMALLOC(xlindx, sizeP1, int);
      CMALLOC(invp, size, int);
      CMALLOC(w3, size, double);

      CMALLOC(cols, spam_zaehler, int);
      CMALLOC(rows, sizeP1, int);
   
      CMALLOC(DD, spam_zaehler, double);
     // prepare spam

     
     
      F77_CALL(spamdnscsr)(&size, &size, M, &size, DD, cols, rows, 
	       &spam_tol); // create spam object
     

      //   for (int k=0; k<spam_zaehler; k++) printf("%f ",M[k]); printf("\n");
      //for (int k=0; k<spam_zaehler; k++) printf("%f ", DD[k]);  printf("\n");
      // for (int k=0; k<spam_zaehler; k++) printf("%d ",cols[k]);  printf("\n");
      //  for (int k=0; k<size+1; k++) printf("%d ",rows[k]);  printf("\n");

     
      // spam solve
      if (rhs_cols <= 0) { // UNBEDINGT VOR double *RHS;
	CMALLOC(RHS, sizeSq, double);	
	RHS_COLS = size;
	for (i=0; i<sizeSq; i += sizeP1) RHS[i] = 1.0;
      } else RHS_COLS = rhs_cols;	
      int totbytes = size * RHS_COLS;
      double *RHS = rhs_cols > 0 ? rhs : pt->RHS;
      

      if (posdef) {
	// calculate spam_cholesky
	err = 4; // to get into the while loop
	while (err == 4 || err == 5) {
	  if (nnzcolindices == 0) {
	    double rel = nnzA / (double) size;
	    if (rel < 5) {
	      nnzcolindices = nnzA * (1.05 * rel - 3.8);
	      if (nnzcolindices < 1000) nnzcolindices = 1000;
	    } else {
	      nnzcolindices = nnzA;
	    }
	    nnzcolindices *= nnzcfact[doperm];
	    if (nnzcolindices < nnzA) nnzcolindices = nnzA;
	  } else if (err == 5) {
	    int tmp = ceil(nnzcolindices * cholincrease_nnzcol);
	    if (PL > PL_IMPORTANT) 
	      PRINTF("Increased 'nnzcolindices' with 'NgPeyton' method\n(currently set to %d from %d)", tmp, nnzR);
	    nnzcolindices = tmp;
	  }
	  if (nnzcolindices < pt->lindx_n) nnzcolindices = pt->lindx_n;

	  if (nnzR == 0) {
	    double u = floor(.4 * pow(nnzA, 1.2));
	    if (u < 4 * nnzA) u = 4 * nnzA;
	    nnzR = u * nnzRfact[doperm];
	  } else if (err == 4) {
	    int tmp = ceil(nnzR * cholincrease_nnzR);
	    if (PL > PL_IMPORTANT) 
	      PRINTF("Increased 'nnzR' with 'NgPeyton' method\n(currently set to %d from %d)", tmp, nnzR);
	    nnzR = tmp;
	  }
	  if (nnzR < pt->lnz_n) nnzR = pt->lnz_n;
	  else if (nnzR > halfsq) nnzR = halfsq;	
	
	  CMALLOC(lindx, nnzcolindices, int);	
	  CMALLOC(lnz, nnzR, double);
	
	  F77_CALL(cholstepwise)(&size, &nnzA, DD, cols, rows, &doperm,
				 invp, pivot, 
				 &nnzlindx, &nnzcolindices, 
				 lindx, // 
				 xlindx,// 
				 &nsuper, // length of lindx
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
	}
   
      
	double *lnz = pt->lnz;
	int *lindx = pt->lindx;
      

	// spam determinant
	if (logdet != NULL) {
	  double tmp = 0.0;
	  for (i=0; i<size; i++) {
	    // printf("%f %d ", lnz[xlnz[i] - 1], xlnz[i] -1);
	    tmp += log(lnz[xlnz[i] - 1]);
	  }
	  //printf(" ok %d \n", size);
	  *logdet = 2.0 * tmp;	  
	}
	

	F77_CALL(backsolves)(&size, &nsuper, &RHS_COLS, lindx, xlindx, 
			     lnz, xlnz, invp, pivot, xsuper, w3, RHS);
		
      }

      if (err != NOERROR || !posdef) {
	// z < - backsolve(a, forwardsolve( t(a),b))
	//	CMALLOC(t_cols, spam_zaehler, int);
	//	CMALLOC(t_rows, sizeP1, int);
	//	CMALLOC(t_DD, spam_zaehler, double);
	//	F77_CALL(transpose)(&size, &size, DD, cols, rows, 
	//			    t_DD. t_cols, t_rows);
	//	F77_CALL(spamforward)(&size, &RHS_COLS, 
	CERR("'spam' failed");
      }

      if (rhs_cols == 0) MEMCOPY(M, RHS, totbytes * sizeof(double));
      if (PL >=  PL_STRUCTURE) PRINTF("'spam' successful\n");
      
      break;
    }
   
    default :
      SERR("unknown method for 'solvePosDef'");
      
    } // switch

    //printf("very end: %d\n", err);
    if (err==NOERROR) break;
  }
 
   //printf("ending: %d\n", err);

   return err; //  -method;
}
  

SEXP solvePosDef(SEXP M, SEXP rhs, SEXP logdet){
  solve_storage pt;
  int err = NOERROR,
    size = ncols(M), 
    rhs_cols = isMatrix(rhs) ? ncols(rhs) : length(rhs)==0 ? 0 : 1;
  SEXP res;

  if (isMatrix(rhs) || rhs_cols==0) {
    PROTECT(res = allocMatrix(REALSXP, size, rhs_cols == 0 ? size : rhs_cols));
  } else {
    PROTECT(res = allocVector(REALSXP, length(rhs)));
  }

  solve_NULL(&pt);
  if (rhs_cols == 0) {
    MEMCOPY(REAL(res), REAL(M), size * size * sizeof(double));
  } else {
    MEMCOPY(REAL(res), REAL(rhs), size * rhs_cols * sizeof(double));
  }
  
  if (!false)
  err = solvePosDef_(REAL(rhs_cols == 0 ? res : M),  size,  true,
		     rhs_cols == 0 ? NULL : REAL(res),  rhs_cols, 
		     length(logdet) == 0 ? NULL : REAL(logdet),
		     &pt, &SolveParam, PL_IMPORTANT);
  // 1e-8: svd_tol = GLOBAL.general.matrixtolerance
  
  //for (int i=0; i< size * rhs_cols; i++) printf("%f ", REAL(res)[i]); printf("B \n");

  solve_DELETE0(&pt);
   // printf("err = %d\n", err);
  // if (err == ERRORM) printf("error: %s\n", ERRORSTRING);
   
  UNPROTECT(1);
  if (err != NOERROR) ERR("'solvePosDef' failed.");
  
  return(res);
}


int invertMatrix(double *M, int size) {
  solve_storage *pt = (solve_storage*) MALLOC(sizeof(solve_storage));
  int err;
  // to do
  err =  solvePosDef_(M, size, false, NULL, 0, NULL, pt, &SolveParam,
		     PL_IMPORTANT);
  solve_DELETE(&pt);
  return err;
}
