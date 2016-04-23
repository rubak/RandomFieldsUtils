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

#include <R_ext/Lapack.h>
#include "RandomFieldsUtils.h"

const char * InversionNames[(int) Diagonal + 1] = {
  "Cholesky",  "SVD",  "SPAM", "QR", "LU", "no method left",
  "method undefined", "Diagonal"};


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
    FREE(x->xja);
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
      x->invp_n = x->cols_n = x->rows_n =x->lindx_n = x->xja_n =
      //
    x->SICH_n = x->MM_n = x->workspaceD_n = x->workspaceU_n = 
    x->VT_n = x->work_n = x->w2_n = x->U_n = x->D_n = x->workLU_n =
      x->lnz_n =  x->DD_n = x->RHS_n = x->w3_n =
    0;

    x->nsuper = x->nnzlindx = x->size = -1;
    x->method = NoInversionMethod;
  
    x->iwork = x->ipiv = 
    x->pivot = x->xlnz = x->snode = x->xsuper = x->xlindx = 
    x->invp = x->cols = x->rows = x->lindx = x->xja =
    NULL;

    x->SICH = x->MM = x->workspaceD = x->workspaceU = 
    x->VT = x->work = x->w2 = x->U = x->D = x->workLU =
    x->lnz = x->DD = x->RHS = x->w3 =
    NULL;
}


int solve2(double *M, int size, bool posdef,
		double *rhs, int rhs_cols,
		double *logdet 		
		){
  assert(size <= 3);
  if (size <= 0) ERR("matrix in 'solvePosDef' of non-positive size.");
  int i;

  double det;
  switch(size){ // Abfrage nach Groesse der Matrix M + Berechnung der Determinante per Hand
  case 1: det = M[0];
    break;
  case 2: det = M[0] * M[3] - M[1] * M[2];
    break;
  case 3: det = 
      M[0] * (M[4] * M[8] - M[5] * M[7]) 
      - M[1] * (M[3] * M[8] - M[5] * M[6]) 
      + M[2] * (M[3] * M[7] - M[4] * M[6]); // Entwicklung nach 1. Spalte
    break;
  default : BUG;
    break;
  }

  if (det == 0 || (posdef && det < 0)) return ERRORFAILED;
  if (logdet != NULL) *logdet = log(det);

  double detinv = 1.0 / det; // determinant of inverse of M
  
  switch(size){
    case 1: // size of matrix == 1
      if (rhs_cols == 0) M[0] = detinv;   
      else for (i=0; i<rhs_cols; rhs[i++] *= detinv);
      break;
    case 2: // size of matrix == 2
      if (rhs_cols == 0) {
	double swap = M[0] * detinv;
	M[0] =  M[3] * detinv;
	M[1] = -M[1] * detinv;
	M[2] = -M[2] * detinv;
	M[3] = swap;
      } else { // rhs_cols != 0
	double *p = rhs,
	  a = M[0] * detinv,
	  b = M[1] * detinv,
	  c = M[2] * detinv,
	  d = M[3] * detinv;
	if (b != 0.0 || c != 0.0) {
	  for (i=0; i<rhs_cols; i++, p+=2) {
	    double swap = d * p[0] - c * p[1];
	    p[1] = a * p[1] - b * p[0];
	    p[0] = swap;
	  }
	} else {
	  for (i=0; i<rhs_cols; i++, p+=2) {
	    double swap = d * p[0];
	    p[1] = a * p[1];
	    p[0] = swap;
	  }
	}
      }
      break;
  case 3: // size of matrix == 3
    if(rhs_cols == 0){ // invert matrix
      double swap0 = detinv * (M[4] * M[8] - M[5] * M[7]),
	swap1 = -detinv * (M[3] * M[8] - M[5] * M[6]),
	swap2 = detinv * (M[3] * M[7] - M[4] * M[6]),
	swap3 = -detinv * (M[1] * M[8] - M[2] * M[7]),
	swap4 = detinv * (M[0] * M[8] - M[2] * M[6]),
	swap5 = -detinv * (M[0] * M[7] - M[1] * M[6]),
	swap6 = detinv * (M[1] * M[5] - M[2] * M[4]),
	swap7 = -detinv * (M[0] * M[5] - M[2] * M[3]);
	M[8] = detinv * (M[0] * M[4] - M[1] * M[3]);
	M[0] = swap0;
	M[1] = swap1;
	M[2] = swap2;
	M[3] = swap3;
	M[4] = swap4;
	M[5] = swap5;
	M[6] = swap6;
	M[7] = swap7;
    } else { // solve system given by M and rhs
      double *p = rhs;
      for (i=0; i<rhs_cols; i++, p+=3) {
	double swap0 = detinv * (   p[0] * (M[4] * M[8] - M[5] * M[7]) 
				  - p[1] * (M[1] * M[8] - M[2] * M[7])
				  + p[2] * (M[1] * M[5] - M[2] * M[4]));
	double swap1 = detinv * (- p[0] * (M[3] * M[8] - M[5] * M[6]) 
				 + p[1] * (M[0] * M[8] - M[2] * M[6]) 
				 - p[2] * (M[0] * M[5] - M[2] * M[3]));
	p[2] = detinv * (   p[0] * (M[3] * M[7] - M[4] * M[6])
			  - p[1] * (M[0] * M[7] - M[1] * M[6]) 
			  + p[2] * (M[0] * M[4] - M[1] * M[3]));
	p[0] = swap0;
	p[1] = swap1;
      }
    }
    break;
  default: BUG;
  }
  
  return NOERROR;
}


int solvePosDef_(double *M, int size, bool posdef,
		 double *rhs, int rhs_cols,
		 double *logdet, 
		 solve_storage *pt, solve_param *Sp, int PL
 		 ){
  /*

    M: a square matrix (symmetry is not checked) of size x size;
       NOTE THAT THE CONTENTS OF M IS DESTROYED IFF NO RHS IS GIVEN
       in case rhs is not given, the inverse of M is returned here
    posdef: whether or not the matrix is positiv (semi)definite --
            to some extend solvePosDef can deal with non-positiv definite
            functions
    rhs : right hand side of the equality with rhs_cols columns
          NOTE THAT THE CONTENTS OF rhs WILL BE DESTROYED IF rhs  GIVEN 
          the solution of the equality is returned in rhs
    logdet : if not NULL the logarithm of the determinant is returned
    pt  : working space. If NULL, internal working spaces are used.
 
          A non-NULL value gives an advantage only if solvePosDef is called
          repeated. Then 
            solve_storage *pt = (solve_storage*) malloc(sizeof(solve_storage);
            solve_NULL(pt);
          prepares pt. Deletion is done only at the very end by
            solve_DELETE(pt);
          In meanwhile the working space is managed by solvePosDef_;
     Sp : parameters. If NULL, the behaviour is as described in
          the R manual for solvePosDef.
          The parameters are described in the package 'RandomFields'
          under ?RFoptions
     PL : printing level; 1 is standard
	  
  */

  // printf("size=%d %f %f %f %f\n", size, M[0], M[1], M[2], M[3]);

  // http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F01/f01intro.xml#
  assert(NA_LOGICAL == INT_MIN && NA_LOGICAL == NA_INTEGER); // nur zur sicherheit, wegen usr_bool
  //          eigentlich sollte usr_bool unabhaengig davon funktionieren

  if (size <= 3) return solve2(M, size, posdef, rhs, rhs_cols, logdet);

  assert(SOLVE_METHODS >= 2);
  solve_param
    sp_standard = solve_param_default,
    *sp = Sp == NULL ? &sp_standard : Sp;
  InversionMethod method = NoInversionMethod;
  int  
    err = NOERROR,
    spam_zaehler = 0,
    nnzA = 0,
    sizeSq = size * size,
    sizeP1 = size + 1;
  usr_bool
    sparse = sp->sparse;
  double 
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
	if (!posdef) for (; j<end; j++) nnzA += fabs(M[j++]) >= spam_tol;
      }
      diag = (nnzA == 0);
      if (posdef) nnzA *= 2;
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
      if (diag && !posdef) {
	for (; j<end; j++) {
	  diag = false;
	  break;
	}
      }
    }
  }
  

  //printf("start vxxx %f\n", (double) sparse);

  if (pt == NULL) BUG;

  if (diag) {
    if (logdet != NULL) {
      double tmp = 0.0;
      for (int i=0; i<sizeSq; i+=sizeP1) tmp += log(M[i]);
      *logdet = tmp;
    }
    if (rhs_cols == 0) {
      for (int i=0; i<sizeSq; i += sizeP1) M[i] = M[i] <= 0.0 ? 0.0 : 1.0/ M[i];
    } else {
      CMALLOC(MM, size, double);
      for (int i=0; i<size; i++) {
	int idx = i * sizeP1;
	MM[i] = M[idx] == 0.0 ? 0.0 : 1.0 / M[idx];
      }
      int j;
      for (int k=j=0; j<rhs_cols; j++)
	for (int i=0; i<size; i++) rhs[k++] *= MM[i];
    }
    err = NOERROR;
    goto ErrorHandling;
  }
  
  //   printf("Ok vxxx %f\n", (double) sparse);

  InversionMethod *Meth;
  if (sparse == True || sp->Methods[0] >= NoFurtherInversionMethod) {
    Meth = pt->newMethods;
    if (sparse == True) {
      pt->newMethods[0] = Sparse;
      pt->newMethods[1] = 
	sp->Methods[0] < NoFurtherInversionMethod && sp->Methods[0] != Sparse
	? sp->Methods[0] 
	: posdef ? Cholesky : LU;  
      if (SOLVE_METHODS > 2) {
	pt->newMethods[2] = 
	  sp->Methods[0]<NoFurtherInversionMethod && sp->Methods[0] != Sparse &&
	  sp->Methods[1]<NoFurtherInversionMethod && sp->Methods[1] != Sparse
	  ? sp->Methods[1] 
	  : posdef ? SVD : LU;
      }
      // pt->newMethods[1] = Sparse;
    } else {
      pt->newMethods[0] = posdef ? Cholesky : LU;  
      pt->newMethods[1] =  posdef ? SVD : LU;
      if (SOLVE_METHODS > 2) pt->newMethods[2] = SVD;
    }
  } else Meth = sp->Methods;

 //  PL = 10;

  if (!posdef && Meth[0] != SVD && Meth[0] != SVD) {
    err = ERRORNOTPROGRAMMEDYET;
    goto ErrorHandling;
  }
 
  if (Meth[1] != Meth[0] && rhs_cols==0) { // at least two different Methods in the list
    CMALLOC(SICH, sizeSq, double);
    memcpy(SICH, M, sizeSq * sizeof(double));
  }
  double *SICH;
  SICH = pt->SICH;

  for (int m=0; m<SOLVE_METHODS && (m==0 || Meth[m] != Meth[m-1]); m++) {
    method = Meth[m];
   if (PL>=PL_STRUCTURE) { 
     PRINTF("method to calculate inverse : %s\n", InversionNames[method]);
    }
 
    // print("m=%d Method %d chol=%d QR=%d SVD=%d sparse=%d\n", m, method, Cholesky, QR, SVD, Sparse);
   
    if (method<0) break;
    if (m > 0 && rhs_cols == 0) {
      memcpy(M, SICH, sizeSq * sizeof(double));
    }

    switch(method) {
 
     case Cholesky : // cholesky       
       if (!posdef) CERR("Cholesky needs positive definite matrix");
       if (size > sp->max_chol) CERR("matrix too large for Cholesky decomposition.");
       if (rhs_cols == 0) {
	 //for (int ii=0; ii<100; ii++) printf("%10.8f ", M[ii]); printf("M\n");
	 F77_CALL(dpotrf)("U", &size, M, &size, &err);  
	 //	   for (int ii=0; ii<100; ii++) printf("%10.8f ", M[ii]); printf("Minv\n");
	 if (err == NOERROR) {
	   int i;
	   if (logdet != NULL) {
	     for (*logdet=0.0, i=0; i < sizeSq; i+=sizeP1) {
	       *logdet += log(M[i]);
	     }
	     *logdet *= 2;
	   }
	   F77_CALL(dpotri)("U", &size, M, &size, &err);  
	   if (rhs_cols == 0) {
	     long  i2, i3, j;
	     for (i2=i=0; i<size; i++, i2+=size + 1) {	
	       for (i3 = i2 + 1, j = i2 + size; j<sizeSq; j+=size) M[i3++]=M[j];
	     }
	   }
	 }
       } else {
	 //printf("sizeSq %d %d\n", sizeSq, pt->MM_n); printf("before\n");
 	 CMALLOC(MM, sizeSq, double);
	 //	 printf("callo done \n\n");
	 MEMCOPY(MM, M, sizeSq * sizeof(double));
	 //	 printf("memco done \n\n");
	 F77_CALL(dpotrf)("U", &size, MM, &size, &err);  
	 if (err == NOERROR) {
	   int i;
	   if (logdet != NULL) {
	     for (*logdet=0.0, i=0; i < sizeSq; i+=sizeP1)
	       *logdet += log(MM[i]);
	     *logdet *= 2;
	   }
 	   F77_CALL(dpotrs)("U", &size, &rhs_cols, MM, &size, rhs, &size, &err);
 	 }
       }
     
       if (err != NOERROR) {	
	 CERR1("Cholesky decomposition failed with err=%d of 'dpotr*' in 'solvePosDef'. Probably matrix not strictly positive definite.\n", err);
       }

       if (PL >=  PL_DETAILSUSER) PRINTF("Cholesky decomposition successful\n");
     
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
      if (PL >=  PL_DETAILSUSER) PRINTF("QR successful\n");
      break;
    }

    case SVD : {// SVD : M = U D VT
      if (size > sp->max_svd) CERR("matrix too large for SVD decomposition.");
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
	int i;
	for (*logdet=0.0, i = 0; i < size; *logdet += log(D[i++]));
      }

      double svd_tol = sp->svd_tol;
      for (int j=0; j<size; j++) D[j] = fabs(D[j]) < svd_tol ? 0.0 : 1.0 / D[j];

      if (rhs_cols > 0) {
	int tot = size * rhs_cols;
	CMALLOC(w2, tot, double);	
	matmulttransposed(U, rhs, w2, size, size, rhs_cols);
	for (k=0; k<tot; )
	  for (int i=0; i<size; i++) w2[k++] *= D[i];
	matmulttransposed(VT, w2, rhs, size, size, rhs_cols);
      } else {
	// calculate inverse of covariance matrix
	int j;
	for (k=0, j=0; j<size; j++) {
	  double dummy = D[j];
	  for (int i=0; i<size; i++) U[k++] *= dummy;
	}
	matmult_tt(U, VT, M, size, size, size); // V * U^T
      } 

      if (PL >=  PL_DETAILSUSER) PRINTF("svd successful\n");
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
 
      if (logdet != NULL) {
	CERR("logdet cannot be determined for 'LU'");
	int i;
	for (*logdet=0.0, i = 0; i < sizeSq; i += sizeP1) *logdet += log(M0[i]);
      }

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
      
      if (PL >=  PL_DETAILSUSER) PRINTF("LU decomposition successful\n");
      break;
    }
     
    case Sparse : {// sparse matrix
     int nnzlindx, RHS_COLS, 	
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
      CMALLOC(w3, size, double);

      CMALLOC(cols, spam_zaehler, int);
      CMALLOC(rows, sizeP1, int);
   
      CMALLOC(DD, spam_zaehler, double);
     // prepare spam

      F77_CALL(spamdnscsr)(&size, &size, M, &size, DD, cols, rows, 
			   &spam_tol); // create spam object
     
      // spam solve
      if (rhs_cols <= 0) { // UNBEDINGT VOR double *RHS;
	CMALLOC(RHS, sizeSq, double); // pt->RHS !!	
	RHS_COLS = size;
	for (int i=0; i<sizeSq; i += sizeP1) RHS[i] = 1.0; // pt->RHS !!
      } else RHS_COLS = rhs_cols;	
      int totbytes = size * RHS_COLS;
      double *RHS = rhs_cols > 0 ? rhs : pt->RHS;
      

      pt->nsuper = 0;
      if (posdef) {
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
				 &nnzlindx, &nnzcolindices, 
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
      
	double *lnz = pt->lnz;
	int *lindx = pt->lindx;
      

	// spam determinant
	if (logdet != NULL) {
	  double tmp = 0.0;
	  for (int i=0; i<size; i++) {
	    tmp += log(lnz[xlnz[i] - 1]);
	  }
	  *logdet = 2.0 * tmp;	  
	}


	F77_CALL(backsolves)(&size, &(pt->nsuper), &RHS_COLS, 
			     lindx, // colindices
			     xlindx, //colpointers
			     lnz, 
			     xlnz, //  rowpointers
			     invp, pivot,
			     xsuper, // supernodes
			     w3, RHS);
		
      } else CERR("'spam' needs a positive definite matrix");

      if (rhs_cols == 0) MEMCOPY(M, RHS, totbytes * sizeof(double));
      if (PL >=  PL_DETAILSUSER) PRINTF("'spam' successful\n");
      
      break;
    }
   
    default :
      GERR("unknown method for 'solvePosDef'");
      
    } // switch

    if (err==NOERROR) break;
  } // for m


 ErrorHandling:
   return err; //  -method;
}
  

SEXP solvePosDef(SEXP M, SEXP rhs, SEXP logdet){
  solve_storage pt;
  int 
    rhs_rows, rhs_cols,
    err = NOERROR,
    size = ncols(M), 
    rows = nrows(M);
  bool deleteMM = false;
  SEXP res;

  if (isMatrix(rhs)) {
    rhs_rows = nrows(rhs);
    rhs_cols = ncols(rhs);
  } else if ((rhs_rows = length(rhs)) == 0) {
    rhs_cols = 0;
  } else {
    rhs_cols = 1;
  }
  if (rows != size) ERR("not a square matrix");
  if (rhs_rows > 0 && rhs_rows != size) ERR("vector size does not match the matrix size");
  
  int 
    new_cols = rhs_cols == 0 ? size : rhs_cols,
    total = size * new_cols;
  if (isMatrix(rhs) || rhs_cols==0) {
    PROTECT(res = allocMatrix(REALSXP, size, new_cols));
  } else {
    PROTECT(res = allocVector(REALSXP, total));
  }

  SEXP from = rhs_cols == 0 ? M : rhs;  
  if (TYPEOF(from) == REALSXP) {
    MEMCOPY(REAL(res), REAL(from), total * sizeof(double));
  } else if (TYPEOF(from) == INTSXP) {
    for (int i=0; i<total; i++) {
      REAL(res)[i] = INTEGER(from)[i] == NA_INTEGER 
	? RF_NA : (double) INTEGER(from)[i];
    }
  } else ERR("numerical matrix and/or vector expected");
  
  double *MM;
  if (rhs_cols == 0) {
    MM = REAL(res);
  } else {
    MM = REAL(M);
    if ((deleteMM = TYPEOF(M) != REALSXP)) {
      if (TYPEOF(from) == INTSXP) {
        MM = (double*) MALLOC(total * sizeof(double));
        for (int i=0; i<total; i++) 
          REAL(M)[i] = INTEGER(M)[i] == NA_INTEGER 
	    ? RF_NA : (double) INTEGER(M)[i];
      }
    }  
  }

  if (size <= 3) {
     err = solve2(MM, size, true, rhs_cols == 0 ? NULL : REAL(res),  rhs_cols, 
      length(logdet) == 0 ? NULL : REAL(logdet));
  } else {
    solve_NULL(&pt);     
    err = solvePosDef_(MM, size, true, rhs_cols == 0 ? NULL : REAL(res), 
                       rhs_cols, length(logdet) == 0 ? NULL : REAL(logdet),
		       &pt, &SolveParam, PL_IMPORTANT);
    solve_DELETE0(&pt);
  }

  if (deleteMM) FREE(MM);
    
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
