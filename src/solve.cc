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

#ifdef _OPENMP
#include <omp.h>
#endif
#include <R_ext/Lapack.h>
#include "RandomFieldsUtils.h"
#include "own.h"

const char * InversionNames[nr_InversionMethods] = {
  "Cholesky",  "SVD",  "SPAM",
  "method undefined",
  "QR", "LU", 
  "no method left", "direct formula", "diagonal"};


    //  double *A_= A, *B_= B;				
     // i_ = N,					

double scalar(double *A, double *B, int N) {
  double ANS;
  SCALAR_PROD(A, B, N, ANS);
  return ANS;
}

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
    FREE(x->w2);
    FREE(x->U);
    FREE(x->D);

    FREE(x->workLU);
  
    FREE(x->lnz); 
    FREE(x->DD);
    FREE(x->w3);
    FREE(x->result);
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
    x->lnz_n =  x->DD_n = x->w3_n = x->result_n =
    0;
  
  x->nsuper = x->nnzlindx = x->size = -1;
  x->method = NoInversionMethod;
  for (int i=0; i<SOLVE_METHODS; x->newMethods[i++] = NoInversionMethod);
  
  x->iwork = x->ipiv = 
    x->pivot = x->xlnz = x->snode = x->xsuper = x->xlindx = 
    x->invp = x->cols = x->rows = x->lindx = x->xja =
    NULL;
  
  x->SICH = x->MM = x->workspaceD = x->workspaceU = 
    x->VT = x->work = x->w2 = x->U = x->D = x->workLU = 
     x->lnz = x->DD = x->w3 = x->result = NULL;
}


int solve3(double *M, int size, bool posdef,
	   double *rhs, int rhs_cols,
	   double *result, double *logdet 		
		){
  assert(size <= 3);
  if (size <= 0) ERR("matrix in 'solvePosDef' of non-positive size.");

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
  case 1 : {// size of matrix == 1
    if (rhs_cols == 0) result[0] = detinv;   
    else for (int i=0; i<rhs_cols; i++) result[i] = rhs[i] * detinv;
  }
    break;
  case 2 : { // size of matrix == 2
    double a = M[0] * detinv,
      d = M[3] * detinv;
    if (rhs_cols == 0) {
      result[0] = d;
      result[1] = -M[1] * detinv;
      result[2] = -M[2] * detinv;
      result[3] = a;
    } else { // rhs_cols != 0
      double *p = rhs, *q = result;
      if (M[1] != 0.0 || M[2] != 0.0) {
	double 
	  b = M[1] * detinv,
	  c = M[2] * detinv;
	for (int i=0; i<rhs_cols; i++, p+=2, q+=2) {
	  double swap = d * p[0] - c * p[1];
	  q[1] = a * p[1] - b * p[0];
	  q[0] = swap;
	}
      } else {
	for (int i=0; i<rhs_cols; i++, p+=2, q+=2) {
	  double swap = d * p[0];
	  q[1] = a * p[1];
	  q[0] = swap;
	}
      }
    }
  }
    break;
  case 3 : {// size of matrix == 3   
    double swap0 = detinv * (M[4] * M[8] - M[5] * M[7]),
      swap1 = detinv * (M[5] * M[6] - M[3] * M[8]),
      swap2 = detinv * (M[3] * M[7] - M[4] * M[6]),
      swap3 = detinv * (M[2] * M[7] - M[1] * M[8]),
      swap4 = detinv * (M[0] * M[8] - M[2] * M[6]),
      swap5 = detinv * (M[1] * M[6] - M[0] * M[7]),
      swap6 = detinv * (M[1] * M[5] - M[2] * M[4]),
      swap7 = detinv * (M[2] * M[3] - M[0] * M[5]),
      swap8 = detinv * (M[0] * M[4] - M[1] * M[3]);
    if(rhs_cols == 0){ // invert matrix
      result[0] = swap0;
      result[1] = swap1;
      result[2] = swap2;
      result[3] = swap3;
      result[4] = swap4;
      result[5] = swap5;
      result[6] = swap6;
      result[7] = swap7;
      result[8] = swap8;
    } else { // solve system given by M and rhs
      double *p = rhs, *q = result;
      for (int i=0; i<rhs_cols; i++, p+=3, q+=3) {
	double swapA = p[0] * swap0 + p[1] * swap3 + p[2] * swap6;
	double swapB = p[0] * swap1 + p[1] * swap4 + p[2] * swap7;
	q[2] = p[0] * swap2 + p[1] * swap5 + p[2] * swap8;
	q[0] = swapA;
	q[1] = swapB;
      }
    }
  }
    break;
  default: BUG;
  }
  
  return NOERROR;
}

int chol3(double *M, int size, double *res){
  // UNBEDINGT in sqrtRHS.cc auch aendern
  assert(size <= 3);
  if (size <= 0) ERR("matrix in 'solvePosDef' of non-positive size.");
  //  if (M[0] < 0) return ERRORFAILED;
  res[0] = sqrt(M[0]);
  if (size == 1) return NOERROR;
  res[1] = 0.0;
  res[size] = M[size] / res[0];
  res[size + 1] = sqrt(M[size + 1] - res[size] * res[size]);
  if (size == 2) return NOERROR;
  res[2] = res[5] = 0.0;
  res[6] = M[6] / res[0];
  res[7] = (M[7] - res[3] * res[6]) / res[4];
  res[8] = sqrt(M[8] - res[6] * res[6] - res[7] * res[7]);
  return NOERROR;
} 


int doPosDef(double *M, int size, bool posdef,
	     double *rhs, int rhs_cols, double *result, double *logdet, 
	     bool sqrtOnly, solve_storage *Pt, solve_param *Sp
	     ){

  /*
    M: (in/out) a square matrix (symmetry is not checked) of size x size;
       NOTE THAT THE CONTENTS OF M IS DESTROYED IFF NO RHS IS GIVEN
       AND result IS NOT GIVEN.
       In case rhs is not given, the inverse of M is returned here
       In case of sqrtonly M is expected to be a positive definite matrix
    posdef (in): whether or not the matrix is positiv (semi)definite --
            to some extend doPosDef can deal with non-positiv definite
            functions
    rhs (in/out) : right hand side of the equality with rhs_cols columns
          NOTE THAT THE CONTENTS OF rhs WILL BE DESTROYED IF rhs IS GIVEN 
 
          the solution of the equality is returned in rhs
    rhs_cols : number of colums of the matrix on the right hand side
    result (out) : NULL or matrix of the size of the result (inverse matrix or 
           of size of the matrix on the right hand side); see also 'M' ans 'rhs'
    logdet (out): if not NULL the logarithm of the determinant is returned
    pt (in/out) : working space. If NULL, internal working spaces are used.
 
          A non-NULL value gives an advantage only if doPosDef is called
          repeatedly. Then 
            solve_storage *pt = (solve_storage*) malloc(sizeof(solve_storage);
            solve_NULL(pt);
          prepares pt. Deletion is done only at the very end by
            solve_DELETE(pt);
          In meanwhile the working space is managed by doPosDef;
     Sp (in): parameters. If NULL, the behaviour is as described in
          the R manual for doPosDef.
          The parameters are described in the package 'RandomFields'
          under ?RFoptions
	  
  */

  
  // http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F01/f01intro.xml#
  assert(NA_LOGICAL == INT_MIN && NA_LOGICAL == NA_INTEGER); // nur zur sicherheit, wegen usr_bool
  //          eigentlich sollte usr_bool unabhaengig davon funktionieren

  double *RESULT = result != NULL ? result : rhs_cols > 0 ? rhs : M;

  //printf("%ld %ld %ld\ %ldn", RESULT, result, rhs, M); BUG;

  if (size <= 3) {
    if (Pt != NULL) {
      Pt->method = direct_formula;
      Pt->size = size;
    }
     return sqrtOnly 
      ? chol3(M, size, RESULT)
      : solve3(M, size, posdef, rhs, rhs_cols, RESULT, logdet);
  }

  assert(SOLVE_METHODS >= 2);
  solve_param
    *sp = Sp == NULL ? &(GLOBAL.solve) : Sp;

  solve_storage *pt;
  if (Pt != NULL) {
    pt = Pt;
    //printf("pt : %ld \n", (long) (pt));    BUG;
 
  } else {
    pt = (solve_storage*) MALLOC(sizeof(solve_storage));
    solve_NULL(pt);    
  }
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
      for (j=i * size; j<end; j++) {
	if (fabs(M[j]) > spam_tol) {
	  diag = false;
	  break;
	}
      }
      if (!diag) break;
      j++;
      end = (i+1) * size;
      if (!posdef) {
	for (; j<end; j++) {
	  diag = false;
	  break;
	}
      }
    }
  }

  
  if (diag) {
    pt->method = Diagonal;
    if (PL>=PL_STRUCTURE) PRINTF("dealing with diagonal matrix\n");
    if (logdet != NULL) {
      double tmp = 0.0;
      for (int i=0; i<sizeSq; i+=sizeP1) tmp += log(M[i]);
      *logdet = tmp;
    }
    if (rhs_cols == 0) {
      MEMCOPY(RESULT, M, sizeSq * sizeof(double));
      if (sqrtOnly) {
	for (int i=0; i<sizeSq; i += sizeP1)
	  RESULT[i] = M[i] > 0.0 ? sqrt(M[i]) : 0.0;	
      } else 
	for (int i=0; i<sizeSq; i += sizeP1) 
	  RESULT[i] = M[i] <= 0.0 ? 0.0 : 1.0 / M[i];
    } else {
      CMALLOC(MM, size, double);
      for (int i=0; i<size; i++) {
	int idx = i * sizeP1;
	MM[i] = M[idx] == 0.0 ? 0.0 : 1.0 / M[idx];
      }
      int j;
      for (int k=j=0; j<rhs_cols; j++)
	for (int i=0; i<size; i++, k++) RESULT[k] = rhs[k] * MM[i];
    }
    err = NOERROR;
    goto ErrorHandling;
  }

  // size of matrix at least 4 x 4, and not diagonal
  InversionMethod *Meth;
  if (sparse == True || sp->Methods[0] == NoFurtherInversionMethod ||
      sp->Methods[0] == NoInversionMethod) {
    Meth = pt->newMethods;
    if (sparse == True) {
      Meth[0] = Sparse;
      bool given0 =  sp->Methods[0] != NoFurtherInversionMethod &&
	sp->Methods[0] != NoInversionMethod;
      Meth[1] = given0 && sp->Methods[0] != Sparse
	? sp->Methods[0] 
	: posdef ? Cholesky : LU;  
      if (SOLVE_METHODS > 2) {
	bool given1 = sp->Methods[1] != NoFurtherInversionMethod &&
	  sp->Methods[1] != NoInversionMethod;
	Meth[2] = given0 && sp->Methods[0] != Sparse &&
	  given1 && sp->Methods[1] != Sparse
	  ? sp->Methods[1] 
	  : posdef ? SVD : LU;
      }
      // pt->newMethods[1] = Sparse;
    } else {
      Meth[0] = posdef ? Cholesky : LU;  
      Meth[1] =  posdef ? SVD : LU;
      if (SOLVE_METHODS > 2) Meth[2] = SVD;
    }
    for (int i=3; i<SOLVE_METHODS; Meth[i++]=NoFurtherInversionMethod);   
  } else Meth = sp->Methods;

  if (!posdef && Meth[0] != SVD && Meth[0] != SVD) {
    err = ERRORNOTPROGRAMMEDYET;
    goto ErrorHandling;
  }

  // cholesky, QR, SVD, LU always destroy original matrix M
  bool gesichert;
  if ((gesichert = rhs_cols==0 && result == NULL)) {
    if ((gesichert = (SOLVE_METHODS > sparse + 1 &&
		      Meth[sparse + 1] != Meth[sparse] &&
		      Meth[sparse + 1] != NoFurtherInversionMethod)
	 || (Meth[sparse] == SVD && sp->svd_tol >= 0.0 && sqrtOnly)
	 )) { // at least two different Methods in the list
      CMALLOC(SICH, sizeSq, double);
      MEMCOPY(SICH, M, sizeSq * sizeof(double));
    }
  }
  double *SICH, *MPT;
  SICH = pt->SICH;

  MPT = M; // also for sparse result
  if (rhs_cols > 0) {
    CMALLOC(MM, sizeSq, double);
    MPT = MM;
  } else if (result != NULL) MPT = result;

 
  //  printf("gesichert %d\n", gesichert);
  //  int size4; size4 = MIN(5, size);printf("MPT\n"); for (int ii=0; ii<size4; ii++) {for (int jj=0; jj<size4; jj++) printf("%e ", M[ii + jj * size]); printf("\n");}; BUG;
 // printf("%ld %ld %ld\n", (long) MPT, (long)M, long(result)); BUG;
 
  //  bool del = GLOBAL.solve.tmp_delete;
  for (int m=0; m<SOLVE_METHODS && (m==0 || Meth[m] != Meth[m-1]); m++) {

    //printf("%d %s\n", m, InversionNames[Meth[m]]); 

    pt->method = Meth[m];
    if (pt->method<0) break;
    if (sqrtOnly) {
      if (pt->method == NoInversionMethod && m<=sparse) BUG;
       if (pt->method == NoFurtherInversionMethod) break;
     if (PL>=PL_STRUCTURE) { 
	PRINTF("method to calculate the square root : %s\n", 
	       InversionNames[pt->method]);
      }
    } else {
	if (PL>=PL_STRUCTURE) { 
	  PRINTF("method to calculate the inverse : %s\n",
		 InversionNames[pt->method]);
      }
    }
     
    if (rhs_cols == 0 && result == NULL) {
      if (m > sparse) {      
	MEMCOPY(MPT, SICH, sizeSq * sizeof(double));
      }
    } else if (pt->method != Sparse) {
      MEMCOPY(MPT, M, sizeSq * sizeof(double));
      
      //           int size2 = MIN(5, size);printf("MPT ueberschrieben\n"); for (int ii=0; ii<size2; ii++) {for (int jj=0; jj<size2; jj++) printf("%e ", MPT[ii + jj * size]); printf("\n");}//BUG;
    }
    
    switch(pt->method) {
    case Cholesky : // cholesky       
      
      if (!posdef) CERR("Cholesky needs positive definite matrix");
      if (size > sp->max_chol)
	CERR("matrix too large for Cholesky decomposition.");
      
      //      if (Sp->cores <= 1 && false) {
      //	F77_CALL(dpotrf)("Upper", &size, MPT, &size, &err);  
      //      } else 
      
      {
	bool VARIABLE_IS_NOT_USED multicore = (GLOBAL.basic.cores > 1);
	//#pragma omp parallel for
        // cmp for instance http://stackoverflow.com/questions/22479258/cholesky-decomposition-with-openmp
	
	err = NOERROR;	  	
	double *A = MPT;
	for (int i=0; i<size; i++, A += size) {
	  double scalar;	  
	  SCALAR_PROD(A, A, i, scalar);
	  if (A[i] <= scalar) {
	    err = ERRORFAILED;
	    break;
	  }
	  A[i] = sqrt(A[i] - scalar);
	  
	  //	  double invsum = 1.0 / A[i];
	  double sum = A[i];
	  double *endB = MPT + sizeSq; 
#ifdef DO_PARALLEL
#pragma omp parallel for if (multicore)
#endif
	  for (double *B=MPT + (i+1) * size; B<endB; B+=size) {
	    double scalar2;
	    SCALAR_PROD(A, B, i, scalar2);
	    //B[i] = invsum * (B[i] - scalar2);
	    B[i] = (B[i] - scalar2) / sum;
	  }
	}
      }


	// see also http://www.wseas.us/e-library/conferences/2013/Dubrovnik/MATHMECH/MATHMECH-25.pdf
	// http://ac.els-cdn.com/0024379586901679/1-s2.0-0024379586901679-main.pdf?_tid=bc5f2c8c-3117-11e6-80e1-00000aab0f02&acdnat=1465789050_ebfe7248d7a126bd2a301e97a3dbf914
	if (false) {
        //braucht 100 % mehr zeit als aufruf von dpotrf
	  // laesst sich nicht ohne weiteres 
	err = NOERROR;	  
	//	  omp_set_num_threads(Sp->cores);
	int isize=0;
	for (int i=0; i<size; i++, isize += size) {
	  double *A = MPT + isize;
	  
	  // #pragma omp parallel for -- warum funktioniert das nicht??
	  for (int j=0; j<i; j++) {
	    int jsize = j * size;
	    double 
	      *B =  MPT + jsize,
	      sum = A[j]; 
	    // printf("a_ij = %f\n", sum);
	    for (int k=0; k<j; k++) sum -= A[k] * B[k];
	    // printf("A_ij = %f %f\n", sum, MPT[jsize + j]);
	    A[j] = sum / B[j];
	  }
	  

	  double sum = A[i] - scalar(A, A, i); 
	  if (sum > 0.0) A[i] = sqrt(sum);
	  else { err = ERRORFAILED; break;}
	}
	}



	if (false) {
        //braucht 100 % mehr zeit als aufruf von dpotrf
	  // laesst sich nicht ohne weiteres 
	err = NOERROR;	  
	//	  omp_set_num_threads(Sp->cores);
	int isize=0;
	for (int i=0; i<size; i++, isize += size) {
	  //#pragma omp parallel for
	  double *A = MPT + isize;
	  //#pragma omp parallel for
	  for (int j=0; j<=i; j++) {
	    int jsize = j * size;
	    double 
	      *B =  MPT + jsize,
	      sum = A[j]; 
	    // printf("a_ij = %f\n", sum);
	    for (int k=0; k<j; k++) sum -= A[k] * B[k];
	    // printf("A_ij = %f %f\n", sum, MPT[jsize + j]);
	    if (j < i) A[j] = sum / B[j];
	    else if (sum > 0.0) A[j] = sqrt(sum);
	    else { err = ERRORFAILED; }

	    // printf("j=%d %f  ", j, MPT[isize + j]);	    
	  }
	}
	}


	if (false) {
	 //https://courses.engr.illinois.edu/cs554/fa2013/notes/07_cholesky.pdf
// saying that no pivoting necessary. needs 150 % more time 
	err = NOERROR;	  
	//	  omp_set_num_threads(Sp->cores);

	for (int k =0; k<size; k++) {
	  int kspalte = k * size; 
	  MPT[k + kspalte] = sqrt(MPT[k + kspalte]);
	  double f = 1.0 / MPT[k + kspalte];
	  for (int i = k + 1; i<size; i++) {
	    int ispalte = i * size;
	    MPT[k + ispalte] *= f;
	  }
#pragma omp parallel for
	  for (int j = k + 1; j<size; j++) {
	    int jspalte = j * size;
	    double factor = MPT[k + jspalte];
	    for (int ispalte=j * size; ispalte<sizeSq; ispalte+=size) {
	      MPT[j + ispalte] -= MPT[k + ispalte] * factor;
	    }
	  }
	}
	}


	
	//	   for (int ii=0; ii<100; ii++) printf("%10.8f ", MPT[ii]); printf("Minv\n");
      if (err == NOERROR) {
	if (sqrtOnly) {	  
	  int deltaend = size - 1;
	  double *end = MPT + sizeSq;
	  for (double *p=MPT + 1; p<end; p+=sizeP1, deltaend--)
	    FILL_IN(p, deltaend, 0.0);
	  
	  

	  /*
	  int deltaend = size;
	  for (int i=0; i<sizeSq; i+=sizeP1) {
	    int end = i + (deltaend--);
	    //	      printf("%d %d\n", i, end);
	    for (int j=i + 1; j<end; MPT[j++]=0.0); // untere Dreiecksmatrix 0
	  }
	  */
	} else {
	  int i;
	  if (logdet != NULL) {
	    for (*logdet=0.0, i=0; i < sizeSq; i+=sizeP1) {
	      *logdet += log(MPT[i]);
	    }
	    *logdet *= 2;
	  }
	  if (rhs_cols == 0) {
	    long  i2, i3, j;
	    F77_CALL(dpotri)("U", &size, MPT, &size, &err);  
	    for (i2=i=0; i<size; i++, i2+=size + 1) {	
	      for (i3 = i2 + 1, j = i2 + size; j<sizeSq; j+=size) 
		MPT[i3++] = MPT[j];
	    }
	  } else {
	    int totalRHS = size * rhs_cols;
	    if (result != NULL) MEMCOPY(RESULT, rhs, sizeof(double) * totalRHS);
	    F77_CALL(dpotrs)("U", &size, &rhs_cols, MPT, &size, RESULT, &size,
			     &err);
	  }
	} // sqrt only
      } // err == NOERROR
      
      if (err != NOERROR) {	
	CERR1("Cholesky decomposition failed with err=%d of 'dpotr*'. Probably matrix not strictly positive definite.\n", err);
      }
      
      if (PL >=  PL_DETAILSUSER) PRINTF("Cholesky decomposition successful\n");
      
      break;
       
    case QR : {// QR returns transposed of the inverse !!
      if (rhs_cols > 0 || logdet != NULL || !sqrtOnly) {
	err = ERRORFAILED;
	continue;
      }

      err = ERRORNOTPROGRAMMEDYET; /// to do: clarify transposed !
      continue;

      CMALLOC(workspaceD, size, double);
      CMALLOC(workspaceU, size, double);

      F77_CALL(dgeqrf)(&size, &size,
		       MPT, &size, // aijmax, &irank, inc, workspaceD, 
		       workspaceU, workspaceD, &size, &err);     
      
      //if (GLOBAL.solve.tmp_delete) {FREEING(workspaceD); FREEING(workspaceU);}
     if (err != NOERROR) {	
	CERR1("'dgeqrf' failed with err=%d\n", err);
      }
      if (PL >=  PL_DETAILSUSER) PRINTF("QR successful\n");
      break;
    }

    case SVD : {// SVD : M = U D VT
      if (size > sp->max_svd) CERR("matrix too large for SVD decomposition.");
      double  optim_lwork, 
        //*Uloc = result,
	*pt_work = &optim_lwork;
      int k = 0,  
	size8 = size * 8,
	lwork = -1;

      CMALLOC(VT, sizeSq, double);
      CMALLOC(U, sizeSq, double);
#define Uloc U
      CMALLOC(D, size, double); 
      CMALLOC(iwork, size8, int);
 
// int size3 = MIN(9, size);printf("MPT\n"); for (int ii=0; ii<size3; ii++) {for (int jj=0; jj<size3; jj++) printf("%e ", MPT[ii + jj * size]); printf("\n");};
 
      for (int i=0; i<=1; i++) {
	F77_CALL(dgesdd)("A", &size, &size, MPT, &size, D, Uloc, &size, VT, &size, 
			 pt_work, &lwork, iwork, &err);
	if (err != NOERROR || ISNAN(D[0])) break;
	lwork = (int) optim_lwork;
	CMALLOC(work, lwork, double);
	pt_work = work;
      }
      if (err != NOERROR) {
	if (PL>PL_ERRORS)
	  PRINTF("Error code F77_CALL(dgesdd) = %d\n", err);
	CERR1("'dgesdd' failed with err=%d\n", err);
	break;
      }
      
//      printf("MPT (%d x %d, %ld %ld M=%ld %ld)=\n", size, size, (long) MPT, (long) Uloc, (long) M, (long) result);int size2 = MIN(5, size);for (int ii=0; ii<size2; ii++) {for (int jj=0; jj<size2; jj++) {printf("%e ", MPT[ii + jj * size]);}printf("\n");}
// printf("Uloc=\n");for (int ii=0; ii<size2; ii++) {for (int jj=0; jj<size2; jj++) {printf("%e ", Uloc[ii + jj * size]);}printf("\n");}printf("\n");
//for (int ii=0; ii<size; ii++) { printf("%e ", D[ii]); }printf("\n");

     if (sqrtOnly) {
	double svdtol = sp->svd_tol;
	/* calculate SQRT of covariance matrix */
	for (int j=0; j<size; j++) {
	  double dummy;
	  dummy = fabs(D[j]) < svdtol ? 0.0 : sqrt(D[j]);
	  for (int i=0; i<size; i++, k++) RESULT[k] = Uloc[k] * dummy;
	}
 
//printf("M=\n");for (int ii=0; ii<size2; ii++) {for (int jj=0; jj<size2; jj++) {printf("%f ", M[ii + jj * size]);}printf("\n");}

// printf("sqrt=\n");for (int ii=0; ii<size2; ii++) {for (int jj=0; jj<size2; jj++) {printf("%f ", RESULT[ii + jj * size]);}printf("\n svdtol=%e\n", svdtol);}

	
	/* check SVD */
 	if (svdtol >= 0.0) {
	  double *Morig = gesichert ? SICH : M;
	  for (int i=0; i<size; i++) {
	    double *Ui = RESULT + i;
	    for (k=i; k<size; k++) {
	       double *Uk = RESULT + k,
		sum = 0.0;
	       for (int j=0; j<sizeSq; j+=size) {
		 sum += Ui[j] * Uk[j];
		 //		 printf("%d %d %d sum=%f %f %f\n", i, k, j, sum, Ui[j], Uk[j]);
	       }
	      
	       //	      printf("i=%d j=%d %f %f\n", i, k, sum, Morig[i * size + k]);
	      
	      if (fabs(Morig[i * size + k] - sum) > svdtol) {
		if (PL > PL_ERRORS) {
		  PRINTF("difference %e at (%d,%d) between the value (%e) of the covariance matrix and the square of its root (%e).\n", 
			 Morig[i * size +k] - sum, i, k, Morig[i*size+k], sum);
		}
		FERR3("required precision not attained  (%e > %e): probably invalid model. See also '%s'\n", fabs(Morig[i * size + k] - sum), svdtol,
		      solve[SOLVE_SVD_TOL]);

		err=ERRORM;
		break;
	      } //else printf("ok (%d,%d) %f %f\n", i, k, Morig[i*size+k],sum);
	    }
	    if (err != NOERROR) break;		
	  }
	  if (err != NOERROR) break;		
	}
  	
      } else {
	// calculate determinant 
	if (logdet != NULL) {
	  double dummy = 0.0;
	  for (int i = 0; i < size; dummy += log(D[i++]));
	  *logdet = dummy;
	}
	
	double svd_tol = sp->svd_tol;
	for (int j=0; j<size; j++) 
	  D[j] = fabs(D[j]) < svd_tol ? 0.0 : 1.0 / D[j];

	if (rhs_cols > 0) {
	  int tot = size * rhs_cols;
	  CMALLOC(w2, tot, double);	
	  matmulttransposed(Uloc, rhs, w2, size, size, rhs_cols);
	  for (k=0; k<tot; )
	    for (int i=0; i<size; i++) w2[k++] *= D[i];
	  matmulttransposed(VT, w2, RESULT, size, size, rhs_cols);
	} else {
	  // calculate inverse of covariance matrix
	  int j;
	  for (k=0, j=0; j<size; j++) {
	    double dummy = D[j];
	    for (int i=0; i<size; i++) Uloc[k++] *= dummy;
	  }
	  matmult_tt(Uloc, VT, RESULT, size, size, size); // V * U^T
	}
      }

      if (PL >=  PL_DETAILSUSER) PRINTF("svd successful\n");
      //      if (GLOBAL.solve.tmp_delete) {FREEING(VT);FREEING(U);FREEING(D);
      //FREEING_INT(iwork);FREEING(work);FREEING(w2);}
      break;
    }

    case LU : {// LU
      if (!sqrtOnly) {
	err = ERRORFAILED;
	continue;
      }
      
      CMALLOC(ipiv, size, int);		    
      F77_CALL(dgetrf)(&size, &size, MPT, &size, ipiv, &err);
      if (err != NOERROR) {
	CERR1("'dgetrf' (LU) failed with err=%d\n", err);
      }
 
      if (logdet != NULL) {
	CERR("logdet cannot be determined for 'LU'");
	int i;
	for (*logdet=0.0, i = 0; i < sizeSq; i += sizeP1) *logdet +=log(MPT[i]);
      }

      if (rhs_cols > 0) {
	int totalRHS = size * rhs_cols;
	if (result != NULL) MEMCOPY(RESULT, rhs, sizeof(double) * totalRHS);
	F77_CALL(dgetrs)("N", &size, &rhs_cols, MPT, &size, ipiv, 
			 RESULT, &size, &err);
	if (err != NOERROR) {	
	  CERR1("'dgetrs' (LU) failed with err=%d\n", err);
	}
      } else {
	int lwork = -1;
	double dummy,
	  *p = &dummy;
	for (int i=0; i<=1; i++) { 
	  F77_CALL(dgetri)(&size, MPT, &size, ipiv, p, &lwork, &err);	
	  if (err != NOERROR) break;
	  lwork = (int) dummy;
	  CMALLOC(workLU, lwork, double);
	  p = workLU;
	}
      }      
      if (PL >=  PL_DETAILSUSER) PRINTF("LU decomposition successful\n");
      //if (GLOBAL.solve.tmp_delete) {FREEING_INT(ipiv);FREEING(workLU);}
      break;
    }
	
    case Sparse : {// sparse matrix
     int nnzlindx, 
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

       if (!posdef) CERR("'spam' needs a positive definite matrix");
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
   
      int nDD = spam_zaehler;
      if (nDD < size) nDD = size;
      CMALLOC(DD, nDD, double);
     // prepare spam

      F77_CALL(spamdnscsr)(&size, &size, M, &size, DD,
			   cols, // ja
			   rows, // ia
			   &spam_tol); // create spam object
     
      pt->nsuper = 0;
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
			       xlnz,  // cols of lnz "ja"
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
      
      // spam solve
      
      if (sqrtOnly) {
	
	//BUG; // unexpected behaviour in spam
	
	nnzR = xlnz[size] - 1;
	CMALLOC(xja, nnzR, int);
	F77_CALL(calcja)(&size, &(pt->nsuper), pt->xsuper, 
			 pt->lindx, pt->xlindx, pt->xlnz, xja);
	for (int i=0; i<size; invp[i++]--); 
	F77_CALL(spamcsrdns)(&size, pt->lnz, xja, pt->xlnz, RESULT);
	for (int i=0; i<size; i++) {
	  int endfor = (i + 1) * size;
	  for (int j = i * (size + 1) + 1; j<endfor; RESULT[j++]=0.0);
	}
	//if (GLOBAL.solve.tmp_delete) {FREEING_INT(pivot); FREEING_INT(xlnz);
	//FREEING_INT(snode); FREEING_INT(xsuper); FREEING_INT(xlindx);
	//FREEING_INT(invp);FREEING(w3);FREEING_INT(cols);FREEING_INT(rows);
	//FREEING(DD);FREEING_INT(lindx);FREEING(lnz);FREEING_INT(xja);}
      } else {     
	double *lnz = pt->lnz;
	int RHS_COLS, 	
	  *lindx = pt->lindx;
	
	// spam determinant
	if (logdet != NULL) {
	  double tmp = 0.0;
	  for (int i=0; i<size; i++) {
	    tmp += log(lnz[xlnz[i] - 1]);
	  }
	  *logdet = 2.0 * tmp;	  
	}
	
	/*   z = .Fortran("backsolves", m = nrow,
	     nsuper, p, a@colindices,
	     a@colpointers, as.double(a@entries),
	     a@rowpointers, a@invpivot, a@pivot,
	     a@supernodes, vector("double",nrow),
	     sol = vector("double",nrow*p),
	     as.vector(b,"double"),
	     NAOK = .Spam$NAOK,PACKAGE = "spam")$sol
	*/
	if (rhs_cols <= 0) { // UNBEDINGT VOR double *RHS;
	  RHS_COLS = size;	
	  FILL_IN(RESULT, sizeSq, 0.0);
	  for (int i=0; i<sizeSq; i += sizeP1) RESULT[i] = 1.0; 
	  
	  // for (int i=0; i<sizeSq; i++) printf("%f ",RESULT[i]); printf("\n");
	  //printf(">> %d %d \n", sizeSq, sizeP1);
	  
	} else {
	  RHS_COLS = rhs_cols;	
	  if (result != NULL) 
	    MEMCOPY(RESULT, rhs, size * rhs_cols * sizeof(double));
	}
	
	//printf("nsuper=%d\n", pt->nsuper);
	//	  for (int ii=0; ii<size; ii++) 
	//printf("%d %d %d %d %e\n", ii, pt->nsuper, sizeP1, xsuper[ii],
	//	   w3[ii]);
	
	//	  if (false)
	//	  for (int jsub=0; jsub<=pt->nsuper; jsub++) {
	//	    int fj = xsuper[1 - 1],
	//	      Lj = xsuper[jsub + 1 - 1] -1;
	//	    printf("%d %d %d\n", jsub, fj, Lj);
	//	    for (int jcol=fj; jcol <= Lj; jcol++) {
	//	      printf("%d,%f  ", jcol, w3[jcol -  1]);
	//	    }	    
	//	  }
	
	//	  for (int jcol=1; jcol <= 600; jcol++) {
	//	    w3[jcol - 1] = jcol;
	//   printf("%d,%f  ", jcol, w3[jcol -  1]);
	//  }	    
	
	
	//	  printf("%ld %ld %d\n", RESULT, rhs, rhs_cols);
	//	  for (int ii=0; ii<size; ii++) printf("%d %e\n", ii, RESULT[ii]);
	//	  BUG;
	
	F77_CALL(backsolves)(&size, &(pt->nsuper), &RHS_COLS, 
			     lindx, // colindices
			     xlindx, //colpointers
			     lnz, 
			     xlnz, //  rowpointers
			     invp, pivot,
			     xsuper, // supernodes
			     w3, RESULT);	
	if (PL >=  PL_DETAILSUSER) PRINTF("'spam' successful\n");
	//if (GLOBAL.solve.tmp_delete) {FREEING_INT(pivot);FREEING_INT(xlnz);
	  //FREEING_INT(snode);FREEING_INT(xsuper);FREEING_INT(xlindx);
	  //FREEING_INT(invp);FREEING(w3);FREEING_INT(cols);FREEING_INT(rows);
	  //FREEING(DD);FREEING_INT(lindx);FREEING(lnz);FREEING_INT(xja);}   
     }
         
      break;
    } // Sparse
   
    default : BUG;
    GERR("unknown method in 'RandomFieldsUtils'");

    } // switch

    if (err==NOERROR) break;
  } // for m


 ErrorHandling:
  if (Pt == NULL) solve_DELETE(&pt);
  //else if (GLOBAL.solve.tmp_delete) {FREEING(SICH); FREEING(MM);}
	  
  return err; //  -method;
}
  
 
SEXP doPosDef(SEXP M, SEXP rhs, SEXP logdet, bool sqrtOnly, 
		 solve_param *Sp){
  int 
    rhs_rows, rhs_cols,
    err = NOERROR,
    size = ncols(M), 
    rows = nrows(M);
  bool deleteMM = false,
    deleteRHS = false;
  SEXP res;
  
  /*
 res = PROTECT(allocMatrix(REALSXP, size, size));
   //  printf("%d %d %d %ld\n", rows, rows*rows, total, (long int) REAL(res)); 
MEMCOPY(REAL(res), REAL(M), size * size * sizeof(double)); F77_CALL(dpotrf)("Upper", &rows, REAL(res), &rows, &err);UNPROTECT(1); return res;
  */


  if (rhs == R_NilValue) {
    rhs_rows = rhs_cols = 0;
  } else if (isMatrix(rhs)) {
    rhs_rows = nrows(rhs);
    rhs_cols = ncols(rhs);
  } else if ((rhs_rows = length(rhs)) == 0) {
    rhs_cols = 0;
  } else {
    rhs_cols = 1;
  }
  if (rows != size) ERR("not a square matrix");
  if (rhs_rows > 0 && rhs_rows != size)
    ERR("vector size does not match the matrix size");
  
  int 
    new_cols = rhs_cols == 0 ? size : rhs_cols,
    total = size * new_cols;

  //  res =  PROTECT(isReal(M) ? duplicate(M): coerceVector(M, REALSXP)); UNPROTECT(1); return res;

  if (rhs_cols==0 || isMatrix(rhs)) {
    res = PROTECT(allocMatrix(REALSXP, size, new_cols));
  } else {
   res =  PROTECT(allocVector(REALSXP, total));
  }


  double *MM=NULL, 
    *RHS = NULL;
  if (TYPEOF(M) != REALSXP) {
    if (TYPEOF(M) != INTSXP && TYPEOF(M) != LGLSXP) 
      GERR("numerical matrix expected");
    if ((deleteMM = rhs_cols != 0))
      MM = (double*) MALLOC(total * sizeof(double));
    else MM = REAL(res);
    if (TYPEOF(M) == INTSXP) {
      for (int i=0; i<total; i++) 
	MM[i] = INTEGER(M)[i] == NA_INTEGER ? RF_NA : (double) INTEGER(M)[i];
    } else {
      for (int i=0; i<total; i++) 
	MM[i] = LOGICAL(M)[i] == NA_LOGICAL ? RF_NA : (double) LOGICAL(M)[i];
    } 
  } else MM = REAL(M); 

  if (rhs_cols > 0) {
    if ((deleteRHS = TYPEOF(rhs) != REALSXP)) {
      if (TYPEOF(res) != INTSXP && TYPEOF(rhs) != LGLSXP) 
	GERR("numerical matrix expected");
      int totalRHS = rhs_cols * rhs_rows; 
      RHS = (double*) MALLOC(totalRHS * sizeof(double));
      if (TYPEOF(rhs) == INTSXP) {
	for (int i=0; i<totalRHS; i++) 
	  RHS[i] = INTEGER(rhs)[i] == NA_INTEGER 
	    ? RF_NA : (double) INTEGER(rhs)[i];
      } else if (TYPEOF(rhs) == LGLSXP) {
	for (int i=0; i<totalRHS; i++) 
	  RHS[i] = LOGICAL(rhs)[i] == NA_LOGICAL
	    ? RF_NA : (double) LOGICAL(rhs)[i];
      } 
    } else RHS = REAL(rhs);
  }
  
  err = doPosDef(MM, size, true, rhs_cols == 0 ? NULL : RHS, rhs_cols, 
		 (rhs_cols == 0 && TYPEOF(M) == REALSXP) ||
		 (rhs_cols > 0 && TYPEOF(rhs) == REALSXP) ? REAL(res) : NULL, 
		 length(logdet) == 0 ? NULL : REAL(logdet),
		 sqrtOnly, NULL, Sp);

 ErrorHandling:
  if (deleteMM) FREE(MM);
  if (deleteRHS) FREE(RHS);
  
  UNPROTECT(1);
  if (err != NOERROR) {
    if (sqrtOnly) {ERR("'cholPosDef' failed");}
    else { ERR("'solvePosDef' failed.");}
  }

  return res;
}


SEXP SolvePosDef(SEXP M, SEXP rhs, SEXP logdet){
  return doPosDef(M, rhs, logdet, false, &(GLOBAL.solve));
}

int solvePosDefResult(double *M, int size, bool posdef, 
		      double *rhs, int rhs_cols, double *result,
		      double *logdet, solve_storage *PT) {
 return doPosDef(M, size, posdef, rhs, rhs_cols, result, logdet, false,
		     PT, &(GLOBAL.solve));
}

int solvePosDef(double *M, int size, bool posdef, 
		 double *rhs, int rhs_cols, 
		 double *logdet, 
		 solve_storage *PT) {
  // printf("solve : %ld %ld %ld %ld\n", (long) (PT), long(M), (long) logdet, (long) rhs); BUG;
  return doPosDef(M, size, posdef, rhs, rhs_cols, NULL, logdet, false,
		  PT, &(GLOBAL.solve));
}


int invertMatrix(double *M, int size) {
  solve_storage *pt = (solve_storage*) MALLOC(sizeof(solve_storage));
  int err;
  // to do
  err =  doPosDef(M, size, false, NULL, 0, NULL, NULL, false,
		  pt, &(GLOBAL.solve));
  solve_DELETE(&pt);
  return err;
}
