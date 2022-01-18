/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2021 Martin Schlather

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


#define RFU_LOCAL 1  


#include "Basic_utils.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include "kleinkram.h"
#include "General_utils.h"
#include "zzz_RandomFieldsUtils.h"
#include "xport_import.h"


const char // constants cannot be exported; 
*KKR_TYPE_NAMES[LAST_R_TYPE_NAME + 1] = { // never change ! see AutoRFU.cc
  "NILSXP" /* 0 */,
  "SYMSXP", "LISTSXP", "CLOSXP", "ENVSXP", "PROMSXP",
  "LANGSXP", "SPECIALSXP", "BUILTINSXP", "CHARSXP", "LGLSXP" /* 10 */,
  "??", "??", "INTSXP", "REALSXP", "CPLXSXP",
  "STRSXP", "DOTSXP", "ANYSXP", "ECSXP", "EXPRSXP" /*20 */,
  "BCODESXP", "EXTPTRSXP", "WEAKREFSXP", "RAWSXP", "S4SXP" /* 25 */,
  "", "", "", "", "NEWSXP" /* 30 */,
  "FREESXP", "??SXP"};



#define USE_OWN_ALG(SCALAR_LEN, PARALLEL) true
#define USE_OWN_SCALAR_PROD true

#define SCALAR(A,B,C) scalarX(A,B,C, SCALAR_AVX)

void strcopyN(char *dest, const char *src, int n) {
  if (n > 1) {
    n--; 
    strncpy(dest, src, n);
  }
  dest[n] = '\0';
}

void AtA(double *a, Long nrow, Long ncol, double *C, int VARIABLE_IS_NOT_USED cores) {
  // C =  A^T %*% A
  if (USE_OWN_ALG(nrow, ncol) || nrow * ncol > MAXINT) {   
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif  
    for (Long i=0; i<ncol; i++) {
      double 
	*A = a + i * nrow,
	*B = A;
      for (Long j=i; j<ncol; j++, B+=nrow) {
      C[i * ncol + j] = C[i + ncol * j] = SCALAR(A, B, nrow);
      }
    }
  } else {
    double alpha = 1.0,
      beta = 0.0;
    int ncol0 = ncol,
      nrow0 = nrow;
    MEMSET(C, 0, ncol * ncol *sizeof(double));
    F77dsyrk("U","T", &ncol0, &nrow0, &alpha, a, &nrow0, &beta, C,
		    &ncol0);
    for (Long i=1; i<ncol; i++) {
      for (Long j=0; j<i; j++) {
	C[i + ncol * j] = C[i * ncol + j];
      }
    }
  }
}
 

void xA_noomp(double *x, double*A, Long nrow, Long ncol, double *y) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y, x, sizeof(double) * nrow);
  } else {
    for (Long i=0; i<ncol; i++) {
      y[i] = SCALAR(x, A + i * nrow, nrow);
    }
  }
}


void xA(double *x, double*A, Long nrow, Long ncol, double *y, int VARIABLE_IS_NOT_USED cores) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y, x, sizeof(double) * nrow);
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (MULTIMINSIZE(ncol) && MULTIMINSIZE(nrow))
#endif  
   for (Long i=0; i<ncol; i++) y[i] = SCALAR(x, A + i * nrow, nrow);
  } 
} 

  
void xA(double *x1, double *x2,  double*A, Long nrow, Long ncol, double *y1,
	double *y2) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y1, x1, sizeof(double) * nrow);
    MEMCOPY(y2, x2, sizeof(double) * nrow);
  } else {
    double *a = A;
    for (Long i=0; i<ncol; i++, a += nrow) {
      y1[i] = SCALAR(x1, a, nrow);
      y2[i] = SCALAR(x2, a, nrow);
    }
  }	
}  
  

double xAx(double *x, double*A, Long nrow, int VARIABLE_IS_NOT_USED cores) {
  if (USE_OWN_SCALAR_PROD || nrow * nrow > MAXINT) {
    double sum = 0.0;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:sum) schedule(static) if (MULTIMINSIZE(nrow) && MULTIMINSIZE(nrow))
#endif  
    for (Long i=0; i<nrow; i++)
      sum += x[i] * SCALAR(x, A + i * nrow, nrow);
    return sum;
  } else {
   double alpha = 1.0,
     beta = 0.0;
    int incx = 1;
    double *y = (double*)  MALLOC(nrow * sizeof(double));
    // z = x^\top A
    int nrow0 = nrow;
    F77dgemv("T", &nrow0, &nrow0, &alpha, A, &nrow0, x, &incx, &beta, y, &incx);
    // z^top x
    alpha = F77ddot(&nrow0, x, &incx, y, &incx);    
    FREE(y);
    return alpha;
  }
}

void Ax(double *A, double*x, Long nrow, Long ncol, double *y, int VARIABLE_IS_NOT_USED cores) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y, x, sizeof(double) * nrow);
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (MULTIMINSIZE(ncol) && MULTIMINSIZE(nrow))
    for (Long j=0; j<nrow; j++) {
      double dummy = 0.0;
      Long k = j;
      for (Long i=0; i<ncol; i++, k+=nrow) { 
	dummy += A[k] * x[i];
      }
      y[j] = dummy;
    }
#else
    for (Long i=0; i<nrow; i++) y[i]=0.0;
    for (Long k=0, i=0; i<ncol; i++) { 
      for (Long j=0; j<nrow; j++) {
	y[j] += A[k++] * x[i];
      }
    }
#endif  
  }
}


void Ax(double *A, double*x1, double*x2, Long nrow, Long ncol, double *y1,
	double *y2) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y1, x1, sizeof(double) * nrow);
    MEMCOPY(y2, x2, sizeof(double) * nrow);
  } else {
    for (Long i=0; i<nrow; i++) y1[i]=y2[i]=0.0;
    for (Long k=0, i=0; i<ncol; i++) { 
      for (Long j=0; j<nrow; j++) {
	y1[j] += A[k] * x1[i];
	y2[j] += A[k++] * x2[i];
      }
    }
  }
}


double XkCXtl(double *X, double *C, Long nrow, Long dim, Long k, Long l,
	      int VARIABLE_IS_NOT_USED cores) {
  // (k-th row of X) * C * (l-th row of X)
  // X is nrow x dim matrix
  // C is dim x dim matrix
  double
    *pX = X + k, 
    *pY = X + l, 
    result = 0.0;
  Long size = nrow * dim;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) reduction(+:result)
#endif
  for (Long j=0; j<size; j+=nrow) {
    double scalar = 0.0;
    Long ci = j * dim;
    for (Long i=0; i<size; i+=nrow) scalar += pX[i] * C[ci++];
    result += scalar * pY[j];
  }
  return result;
}


void XCXt(double *X, double *C, double *V, Long nrow, Long dim /* dim of C */, int VARIABLE_IS_NOT_USED cores) {
  Long size = nrow * dim;
  double  
    *endpX = X + nrow,
    *dummy = (double*) MALLOC(sizeof(double) * size); // dummy = XC
  if (dummy == NULL) RFERROR("XCXt: memory allocation error in XCXt");
 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores)  schedule(static)
#endif
  for (double *pX = X; pX < endpX; pX++) {
    double *pdummy = dummy + (pX - X);
    for (Long ci=0, cd=0; cd<size; cd+=nrow) {
      double scalar = 0.0;
      for (Long i=0; i<size; i+=nrow) {
        scalar += pX[i] * C[ci++];
      }
      pdummy[cd] = scalar;
    }
  }

  // V = dummy X^t
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores)  schedule(dynamic, 20)
#endif
  for (Long rv=0; rv<nrow; rv++) {
    for (Long cv=rv; cv<nrow; cv++) {
      double scalar=0.0;
      for (Long i=0; i<size; i+=nrow) {
	scalar += dummy[rv + i] * X[cv + i];
     }
      V[rv + cv * nrow] = V[cv + rv * nrow] = scalar;
    }
  }

  UNCONDFREE(dummy);
}


double xUy(double *x, double *U, double *y, Long dim, int VARIABLE_IS_NOT_USED cores) {
  // U a symmetric matrix given by its upper triangular part
  double xVy = 0.0;
  Long    dimM1 = dim - 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (MULTIMINSIZE(dim)) reduction(+:xVy) 
#endif  
  for (Long d=0; d<dim; d++) {
    Long i, 
      j = dim * d;
    double dummy = 0.0;
    for (i=0; i<=d; i++) dummy += x[i] * U[j++];
    for (j += dimM1; i<dim; i++, j+=dim) dummy += x[i] * U[j];
    xVy += dummy * y[d];
  }
  return xVy;
}

/*

  // U a symmetric matrix given by its upper triangular part
  assert(z != NULL);
  Long   dimM1 = dim - 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (MULTIMINSIZE(dim))
#endif  
  for (Long d=0; d<dim; d++) {
    double dummy;
    Long i,
      j = dim * d;
    for (dummy = 0.0, i=0; i<=d; i++) dummy += x[i] * U[j++];
    for (j += dimM1; i<dim; i++, j+=dim) dummy += x[i] * U[j];
    if (z!=NULL) z[d] = dummy;
  }
  double xVx;
  SCALAR_PROD(z, x, dim, xVx);
  return xVx;

 */

double xUxz(double *x, double *U, Long dim, double *z, int VARIABLE_IS_NOT_USED cores) {
 double xVx = 0.0;
  Long dimM1 = dim - 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) reduction(+:xVx)
#endif
  for (Long d=0; d<dim; d++) {
    Long i, 
      j = dim * d;
    double dummy = 0.0;
    for (dummy = 0.0, i=0; i<=d; i++) dummy += x[i] * U[j++];
    for (j += dimM1; i<dim; i++, j+=dim) dummy += x[i] * U[j];
    if (z != NULL) z[d] = dummy;
    xVx += dummy * x[d];
  }
  return xVx;
}

double xUx(double *x, double *U, Long dim, int VARIABLE_IS_NOT_USED cores) {
  return xUxz(x, U, dim, NULL, cores);
}

double x_UxPz(double *x, double *U, double *z, Long dim, int VARIABLE_IS_NOT_USED cores) {
// x^t (Ux + z); U dreieckmatrix
  double xVx = 0.0;
  Long    dimM1 = dim - 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) reduction(+:xVx)
#endif
  for (Long d=0; d<dim; d++) {
    Long i,
      j = dim * d;
    double dummy = z[d];
    for (i=0; i<=d; i++) dummy += x[i] * U[j++];
    for (j += dimM1; i<dim; i++, j+=dim) dummy += x[i] * U[j];
    xVx += dummy * x[d];
  }
  return xVx;
}



void matmult(double *a, double *b, double *c, Long l, Long m, Long n,
	     int VARIABLE_IS_NOT_USED cores) {
// multiplying an lxm- and an mxn-matrix, saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
   for (Long i=0; i<l; i++) {
     double *A = a + i,
       *C = c + i;
     for (Long j=0; j<n; j++) {
       double dummy = 0.0,
	 *B = b + j * m;
       for (Long k=0; k<m; k++) dummy += A[k*l] * B[k];
       C[j * l] = dummy;
     }
   }
}


double *matrixmult(double *m1, double *m2, Long dim1, Long dim2, Long dim3,
		   int VARIABLE_IS_NOT_USED cores) {
  double *m0 = (double*) MALLOC(sizeof(double) * dim1 * dim3);
  matmult(m1, m2, m0, dim1, dim2, dim3, cores);
  return m0;
}


void Xmatmult(double *A, double *B, double *C, Long l, Long m, Long n,
	      int VARIABLE_IS_NOT_USED cores) {
// multiplying an lxm- and an mxn-matrix, saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<l; i++) {
    for (Long jl=i, jm=0, j=0; j<n; j++, jl+=l, jm+=m) {
      double dummy = 0.0;
      Long endfor = jm + m;
      for (Long kl=i, k=jm; k<endfor; k++, kl+=l) dummy += A[kl] * B[k]; 
      C[jl] = dummy;
    }
  }
}

void matmulttransposed(double *A, double *B, double *c, Long m, Long l, Long n,
		       int VARIABLE_IS_NOT_USED cores) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<l; i++) {    
    double *C = c + i,
      *Aim = A + i * m;
    for (Long j=0; j<n; j++) C[j * l] = SCALAR(Aim, B + j * m, m);
  }
}


/*
void matmulttransposedInt(int *A, int *B, int *c, Long m, Long l, Long n) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<l; i++) {    
    int *C = c + i,
      *Aim = A + i * m;
    for (Long j=0; j<n; j++) C[j * l] = SCALARINT(Aim, B + j * m, m);
  }
}

*/



void matmult_2ndtransp(double *a, double *B, double *c, Long l, Long m, Long n,
		       int VARIABLE_IS_NOT_USED cores) {
// multiplying A and t(B) with dim(A)=(l, m) and dim(B)=(n, m),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (l * m * n > 1000)
#endif
  for (Long i=0; i<l; i++) {
    double *C = c + i,
      *A = a + i;
    for (Long j=0; j<n; j++) {
       double dummy = 0.0,
	 *Bj = B + j;
       for (Long k=0; k<m; k++) dummy += A[k * l] * Bj[k * n];
       C[j*l] = dummy;
    }
  }
}


void matmult_2ndtransp(double *a, double *B, double *c, Long l, Long m,
		       int VARIABLE_IS_NOT_USED cores) {
// multiplying A and t(B) with dim(A)=(l, m) and dim(B)=(l, m),
// saving result in C
  Long lm = l  * m;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (l * m * l > 1000)
#endif
  for (Long i=0; i<l; i++) {
    double *C = c + i,
      *A = a + i;
    for (Long j=0; j<l; j++) {
       double dummy = 0.0,
	 *Bj = B + j;
       for (Long k=0; k<lm; k+=l) dummy += A[k] * Bj[k];
       C[j*l] = dummy;
    }
  }
}



void Xmatmulttransposed(double *A, double *B, double *C, Long m, Long l, Long n,
			int VARIABLE_IS_NOT_USED cores) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<l; i++) {
    Long im = i * m;
    for (Long jl=i, jm=0, j=0; j<n; j++, jl+=l, jm+=m) {
      double dummy = 0.0;
      Long endfor = im + m;
      for (Long jmk=jm, k=im; k<endfor; k++) dummy += A[k] * B[jmk++]; 
      C[jl] = dummy;
    }
  }
}



void matmult_tt(double *a, double *B, double *c, Long m, Long l, Long n,
		int VARIABLE_IS_NOT_USED cores) {
// calculating t(A B) with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<l; i++) {
    double *A = a + i,
      *C = c + i * l;
    for (Long j=0; j<n; j++) {
      double dummy = 0.0,
	*Bjm = B + j * m;
      for (Long k=0; k<m; k++) dummy += A[k * l] * Bjm[k];
      C[j] = dummy;
    }
  }
}



SEXP TooLarge(int *n, Long l){
#define nTooLarge 2 // mit op
  const char *tooLarge[nTooLarge] = {"size", "msg"};
  SEXP namevec, msg;
  PROTECT(msg=allocVector(VECSXP, nTooLarge));
  PROTECT(namevec = allocVector(STRSXP, nTooLarge));
  for (Long i=0; i<nTooLarge; i++)
    SET_STRING_ELT(namevec, i, mkChar(tooLarge[i]));
  setAttrib(msg, R_NamesSymbol, namevec);
  Long i=0;
  SET_VECTOR_ELT(msg, i++, Int(n, l));
  SET_VECTOR_ELT(msg, i,
		 mkString("too many elements - increase max.elements"));
  UNPROTECT(2);
  return msg;
}
SEXP TooLarge(Long n){int nn=(int) n;  return TooLarge(&nn, 1); }
SEXP TooLarge(Long row, Long col){
   int nn[2] = {(int) row, (int) col};
   return TooLarge(nn, 2);
}
 
SEXP TooSmall(){
  SEXP namevec;
  const char *msg = "value has not been initialized";
  PROTECT(namevec = allocVector(STRSXP, 1));
  SET_STRING_ELT(namevec, 0, mkChar(msg));
  UNPROTECT(1);
  return namevec;
}


SEXP Int(int *V, Long n, Long max) {
  SEXP dummy;
  if (V==NULL) return allocVector(INTSXP, 0);
  if (n>max) return TooLarge(n);
   if (n<0) return TooSmall();
   PROTECT(dummy=allocVector(INTSXP, n));
  for (Long i=0; i<n; i++) INTEGER(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Int(int* V, Long n) { return Int(V, n, n); }


SEXP Logic(bool* V, Long n, Long max) {
  SEXP dummy;
  if (V==NULL) return allocVector(VECSXP, 0);
  if (n>max) return TooLarge(n);
  if (n<0) return TooSmall();
  PROTECT(dummy=allocVector(LGLSXP, n));
  for (Long i=0; i<n; i++) LOGICAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}
SEXP Logic(bool* V, Long n) { return Logic(V, n, n); }

SEXP Num(double* V, Long n, Long max) {
  SEXP dummy;
  if (V==NULL) return allocVector(REALSXP, 0);
  if (n>max) return TooLarge(n);
   if (n<0) return TooSmall();
  PROTECT(dummy=allocVector(REALSXP, n));
  for (Long i=0; i<n; i++) REAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}
SEXP Num(double* V, Long n) {  return Num(V, n, n); }

SEXP Result(double* V, Long n, Long max) {
  SEXP dummy;
  if (V==NULL) return allocVector(REALSXP, 0);
  if (n>max) return TooLarge(n);
   if (n<0) return TooSmall();
 PROTECT(dummy=allocVector(REALSXP, n));
  for (Long i=0; i<n; i++) REAL(dummy)[i] = (double) V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Result(double* V, Long n) {
  return Result(V, n, n);
}

SEXP Char(const char **V, Long n, Long max) {
  SEXP dummy;
  if (V==NULL) return allocVector(STRSXP, 0);
  if (n>max) return TooLarge(n);
   if (n<0) return TooSmall();
   PROTECT(dummy=allocVector(STRSXP, n));
   for (Long i=0; i<n; i++){
     SET_STRING_ELT(dummy, i, mkChar(V[i]));  
   }
  UNPROTECT(1);
  return dummy;
}

SEXP Char(const char **V, Long n) { return Char(V, n, n); }

SEXP Mat(double* V, Long row, Long col, Long max) {
  if (V==NULL) return allocMatrix(REALSXP, 0, 0);
  Long n = row * col;
  if (n>max) return TooLarge(row, col);
  SEXP dummy;
  PROTECT(dummy=allocMatrix(REALSXP, row, col));
  for (Long i=0; i<n; i++) REAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Mat(double* V, Long row, Long col) {
  return Mat(V, row, col, MAXINT);
}


SEXP Mat_t(double* V, Long row, Long col, Long max) {
  if (V==NULL) return allocMatrix(REALSXP, 0, 0);
  Long n = row * col;
  if (n>max) return TooLarge(row, col);
  SEXP dummy;
  PROTECT(dummy=allocMatrix(REALSXP, row, col));
  for (Long k=0, j=0; j<col; j++) {
     for (Long i=0; i<row; i++) {
      REAL(dummy)[k++] = V[j + col * i];
    }
  }
  UNPROTECT(1);
  return dummy;
}

SEXP Mat_t(double* V, Long row, Long col) {
  return Mat_t(V, row, col, MAXINT);
}


SEXP MatString(char **V, Long row, Long col, Long max) {
  if (V==NULL) return allocMatrix(STRSXP, 0, 0);
  Long n = row * col;
  if (n>max) return TooLarge(row, col);
  SEXP dummy;
  PROTECT(dummy=allocMatrix(STRSXP, row, col));
  for (Long k=0; k<n; k++)
    SET_STRING_ELT(dummy, k, mkChar(V[k]));
  UNPROTECT(1);
  return dummy;
}

SEXP MatString(char** V, Long row, Long col) {
  return MatString(V, row, col, MAXINT);
}

SEXP MatInt(int* V, Long row, Long col, Long max) {
  if (V==NULL) return allocMatrix(INTSXP, 0, 0);
  Long n = row * col;
  if (n>max) return TooLarge(row, col);
  SEXP dummy;
  PROTECT(dummy=allocMatrix(INTSXP, row, col));
  for (Long i=0; i<n; i++) INTEGER(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP MatInt(int* V, Long row, Long col) {
  return MatInt(V, row, col, MAXINT);
}

SEXP Array3D(double** V, Long depth, Long row, Long col, Long max) {
  if (V==NULL) return alloc3DArray(REALSXP, 0, 0, 0);
  Long
    m = row * col,
    n = row * col * depth;
  if (n>max) {
    int nn[3] = { (int) row, (int) col, (int) depth };
    return TooLarge(nn, 3);
  }
  SEXP dummy;
  PROTECT(dummy=alloc3DArray(REALSXP, depth, row, col));
  for (Long j=0; j<depth; j++) {
    for (Long i=0; i<m; i++) {
      REAL(dummy)[j*m+i] = V[j][i];
    }
  }
  UNPROTECT(1);
  return dummy;
}

SEXP Array3D(double** V, Long depth, Long row, Long col) {
  return Array3D(V, depth, row, col, MAXINT);
}




usr_bool UsrBoolRelaxed(SEXP p, char *name, Long idx) {
  double dummy = Real(p, name, idx);
  if (!R_finite(dummy)) return Nan;
  return dummy==0.0 ? False : True ;
}


usr_bool UsrBool(SEXP p, char *name, Long idx) {
  double dummy = Real(p, name, idx);
  if (dummy == 0.0) return False;
  else if (dummy == 1.0) return True;
  else if (ISNAN(dummy)) return Nan;
  RFERROR2("invalid value (%d) for boolean variable '%.50s'.", (int) dummy, name);
}




SEXP String(char *V) {
  SEXP str;
  PROTECT(str = allocVector(STRSXP, 1)); 
  SET_STRING_ELT(str, 1, mkChar(V));
  UNPROTECT(1);
  return str;
}

SEXP String(char V[][MAXCHAR], Long n, Long max) {
  SEXP str;
  if (V==NULL) return allocVector(STRSXP, 0);
  if (n>max) return TooLarge(n);
  if (n<0) return TooSmall();
  PROTECT(str = allocVector(STRSXP, n)); 
  for (Long i=0; i<n; i++) {
    SET_STRING_ELT(str, i, mkChar(V[i]));
  }
  UNPROTECT(1);
  return str;
}

SEXP String(char V[][MAXCHAR], Long n) {return String(V, n, n);}
 

SEXP String(int *V, const char * List[], Long n, Long endvalue) {
  SEXP str;
  if (V==NULL || n <= 0) return allocVector(STRSXP, 0);
  Long k;
  for (k=0; k<n; k++) {
    //    printf("k=%d %d: %d %d\n", k,n, V[k], NoInversionMethod);
    assert(V[k] <= endvalue);
    if (V[k] == endvalue) break;
  }
  PROTECT(str = allocVector(STRSXP, k)); 
  for (Long i=0; i<k; i++) {
    //printf("i=%d %d %d\n", i,k, V[i]);
    //    printf("i=%d %d %s\n", i,k, List[V[i]]);
    SET_STRING_ELT(str, i, mkChar(List[V[i]]));
  }
  //  printf("ense\n");
  UNPROTECT(1);
  return str;
}

double Real(SEXP p, char *name, Long idx) {
  if (p != R_NilValue) {
    assert(idx < length(p));
    switch (TYPEOF(p)) {
    case REALSXP :  return REAL(p)[idx];
    case INTSXP :
      if (INTEGER(p)[idx]==NA_INTEGER) return RF_NA;
      else return((double) INTEGER(p)[idx]);
    case LGLSXP :
      if (LOGICAL(p)[idx]==NA_LOGICAL) return(RF_NA);
      else return((double) LOGICAL(p)[idx]);
    default : {}
    }
  }

  RFERROR2("'%.50s' can not be transformed to double! (got'%.50s')\n",
	   name,
	   TYPEOF(p) <= LAST_R_TYPE_NAME ? KKR_TYPE_NAMES[TYPEOF(p)] :
	   "something else"
	   );  
  return RF_NA;  // to avoid warning from compiler
}



void Real(SEXP el,  char *name, double *vec, Long maxn) {
  if (el == R_NilValue) {
    RFERROR1("'%.50s' cannot be transformed to double.\n", name);
  }
  Long n = length(el);
  for (Long j=0, i=0; i<maxn; i++) {
    vec[i] = Real(el, name, j);
    if (++j >= n) j=0;
  }
  return;
}


int Integer(SEXP p, char *name, Long idx, bool nulltoNA) {
  //printf("integer %s %d %d len=%d\n",  name, idx, nulltoNA, length(p));
  if (p != R_NilValue) {
    assert(idx < length(p));
    switch(TYPEOF(p)) {
    case INTSXP : 
      return INTEGER(p)[idx]; 
    case REALSXP : 
      double value;
      value = REAL(p)[idx];
      if (ISNAN(value)) {
	return NA_INTEGER;
      }
      int intvalue;
      intvalue = (int) value;
      if (value == intvalue) return intvalue;      
      else {
	RFERROR2("%.50s: integer value expected. Got %10e.", name, value);
      }
    case LGLSXP :
      if (LOGICAL(p)[idx]==NA_LOGICAL) return(NA_INTEGER);
      else return((int) LOGICAL(p)[idx]);
    default : {}
    }
  } else if (nulltoNA) return NA_INTEGER;

  RFERROR2("%.50s: incorrect type. Got '%.50s'.",
	   name,
	   TYPEOF(p) <= LAST_R_TYPE_NAME ? KKR_TYPE_NAMES[TYPEOF(p)]
	   : "something else");

 return NA_INTEGER; // compiler warning vermeiden
}

int Integer(SEXP p, char *name, Long idx) {
  return Integer(p, name, idx, false);
}


void Integer(SEXP el, char *name, int *vec, Long maxn) {
  if (el == R_NilValue) {
    RFERROR1("'%.50s' cannot be transformed to integer.\n",name);
  }
  Long n = length(el);
  for (Long j=0, i=0; i<maxn; i++) {
    vec[i] = Integer(el, name, j);
    if (++j >= n) j=0;
  }
}




void Integer2(SEXP el, char *name, int *vec) {
  Long n;
  if (el == R_NilValue || (n = length(el))==0) {
      RFERROR1("'%.50s' cannot be transformed to integer.\n",name);
  }
 
  vec[0] = Integer(el, name, 0);
  if (vec[0] != NA_INTEGER && vec[0] < 1)
    RFERROR1("first component of '%.50s' must be at least 1", name);
  if (n==1) vec[1] = vec[0];
  else {
    vec[1] = Integer(el, name, n-1);    
    if ( vec[1] != NA_INTEGER && vec[1] < vec[0])
      RFERROR1("'%.50s' must be increasing", name);
    if (n > 2) {
      vec[1] = vec[0];
      for (Long i = 1; i<n; i++)
	if (Integer(el, name, i) != ++(vec[1]))
	  RFERROR1("'%.50s' is not a sequence of numbers",name);
    }
  }
}





bool Logical(SEXP p, char *name, Long idx) {
   if (p != R_NilValue)
    assert(idx < length(p));
   switch (TYPEOF(p)) {
    case REALSXP:
      if (ISNAN(REAL(p)[idx])) return NA_LOGICAL ;
      else return (bool) REAL(p)[idx];
    case INTSXP :
      if (INTEGER(p)[idx]==NA_INTEGER) return NA_LOGICAL;
      else return (bool) INTEGER(p)[idx];
    case LGLSXP : return LOGICAL(p)[idx];
    default : {}
    }
  RFERROR1("'%.50s' cannot be transformed to logical.\n", name);  
  return NA_LOGICAL;  // to avoid warning from compiler
}


char Char(SEXP el, char *name) {
  SEXPTYPE type;
  if (el == R_NilValue) goto ErrorHandling;
  type = TYPEOF(el);
  if (type == CHARSXP) return CHAR(el)[0];
  if (type == STRSXP) {
    if (length(el)==1) {
      if (STRLEN(CHAR(STRING_ELT(el,0))) == 1)
	return (CHAR(STRING_ELT(el,0)))[0];
      else if (STRLEN(CHAR(STRING_ELT(el,0))) == 0)
	return '\0';
    }
  }
 
 ErrorHandling:
  RFERROR1("'%.50s' cannot be transformed to character.\n",  name);  
  return 0; // to avoid warning from compiler
}


void String(SEXP el, char *name, char names[][MAXCHAR], Long maxlen) {
  Long l = length(el);
  SEXPTYPE type;  
  if (el == R_NilValue) goto ErrorHandling;
  if (l > maxlen)  {
    RFERROR1("number of variable names exceeds %d. Take abbreviations?",
	     (int) maxlen);
  }
  type = TYPEOF(el);
  if (type == CHARSXP) {
    for (Long i=0; i<l; i++) {
      names[i][0] = CHAR(el)[i];
      names[i][1] = '\0';
    }
  } else if (type == STRSXP) {
    for (Long i=0; i<l; i++) {
      strcopyN(names[i], CHAR(STRING_ELT(el, i)), MAXCHAR);
    }
  } else goto ErrorHandling;
  return;
 
 ErrorHandling:
  RFERROR1("'%.50s' cannot be transformed to character.\n",  name);  
}


int NonNegInteger(SEXP el, char *name) {
  int num = INT;
  if (num<0) {
    num=0; 
    WARN1("'%.50s', which has been negative, is set 0.\n",name);
  }
  return num;
}

double NonNegReal(SEXP el, char *name) {
  double num = NUM;
  if (num<0.0) {
    num=0.0; 
    WARN1("%.50s, which has been negative, is set 0.\n",name);
   }
  return num;
}

double NonPosReal(SEXP el, char *name) {
  double num = NUM;
  if (num>0.0) {
    num=0.0; 
    WARN1("%.50s, which has been positive, is set 0.\n",name);
  }
  return num;
}

int PositiveInteger(SEXP el, char *name) {
  int num = INT;
  if (num<=0) {
    WARN2("'%.50s', which has been %.50s, is set 1.\n",
	  name, num ? "negative" : "0");
    num=1;
  }
   return num;
}

double PositiveReal(SEXP el, char *name) {
  double num = NUM;
  if (num<=0.0) {
     WARN2("'%.50s', which has been %.50s, is set 1.\n",
	   name, num==0.0 ? "0" : "negative");
     num=1.0; 
   }
  return num;
}



SEXP ExtendedInteger(double x) {
  return ScalarInteger(R_FINITE(x) ? x : NA_INTEGER);
}

SEXP ExtendedBooleanUsr(usr_bool x) {
  return ScalarLogical((int) x);
}


int Match(char *name, name_type List, int n) {
  // == NOMATCHING, -1, if no matching function is found
  // == MULTIPLEMATCHING,-2, if multiple matching fctns are found,  
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  unsigned int ln;
  int Nr;
  Nr=0;
  ln=STRLEN(name);

  while ( Nr < n  && STRNCMP(name, List[Nr], ln)) {
    Nr++;
  }
  if (Nr < n) { 
    if (ln==STRLEN(List[Nr])) // exactmatching -- take first -- changed 1/7/07
      return Nr;
    // a matching function is found. Are there other functions that match?
    int j; 
    bool multiplematching=false;
    j=Nr+1; // if two or more covariance functions have the same name 
    //            the last one is taken 
    while (j<n) {
      while ( (j<n) && STRNCMP(name, List[j], ln)) {j++;}
      if (j<n) {
	if (ln==STRLEN(List[j])) { // exactmatching -- take first 
	  return j;
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) {return MULTIPLEMATCHING;}
  } else return NOMATCHING;
  return Nr;
}

int Match(char *name, const char * List[], int n) {
  // printf("Matching\n");
   // == -1 if no matching name is found
  // == -2 if multiple matching names are found, without one matching exactly
  unsigned int ln;
  int Nr;
  Nr=0;
  ln=STRLEN(name);

  while ( Nr < n  && STRNCMP(name, List[Nr], ln)) {
    //   printf("%.50s\n", List[Nr]);
    Nr++;
  }
  if (Nr < n) { 
    if (ln==STRLEN(List[Nr])) {// exactmatching -- take first -- changed 1/7/07
      return Nr;
    }
    // a matching function is found. Are there other functions that match?
    int j; 
    bool multiplematching=false;
    j=Nr+1; // if two or more covariance functions have the same name 
    //            the last one is taken 
    while (j<n) {
      while ( (j<n) && STRNCMP(name, List[j], ln)) {j++;}
      if (j<n) {
	if (ln==STRLEN(List[j])) { // exactmatching -- take first 
	  return j;
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) {return MULTIPLEMATCHING;}
  } else return NOMATCHING;

  return Nr;
}



void GetName(SEXP el, char *name, const char * List[], int n,
	     int defaultvalue, int endvalue, int *ans, int maxlen_ans) {
  char dummy[1000];
  int
    k = 0, // globale variable !
    len_el = length(el);

  if (len_el > maxlen_ans) 
    RFERROR2("option '%.50s' is too lengthy. Maximum length is %d.",
	     name, maxlen_ans);
  
  if (TYPEOF(el) == STRSXP) {    
    for ( ; k<len_el; k++) {
      ans[k] = Match((char*) CHAR(STRING_ELT(el, k)), List, n);
      if (ans[k] < 0) {
	if (STRCMP((char*) CHAR(STRING_ELT(el, k)), " ") == 0 ||
	    STRCMP((char*) CHAR(STRING_ELT(el, k)), "") == 0) {
	  goto ErrorHandling;
	}
	goto ErrorHandling0;
      }
    }
  } else {
    Integer(el, name, ans, maxlen_ans);
    for (k=0; k<len_el; k++)
      if (ans[k] < 0 || ans[k] >= n) goto ErrorHandling0;
  }
 
  for (k=len_el; k<maxlen_ans; k++) ans[k] = endvalue;
  return;

ErrorHandling0:
  if (TYPEOF(el) == STRSXP)
    SPRINTF(dummy, "'%.50s': unknown value '%.50s'. Possible values are:", 
	    name, CHAR(STRING_ELT(el, k)));
  else 
     SPRINTF(dummy,
	     "'%.50s':  value '%d' not in {0,...%d}. Other possible values are:", 
	     name, ans[k], n-1);
  
  //  printf("'%s': ", CHAR(STRING_ELT(el, k)));
  int i;
  for (i=0; i<n-1; i++) {
    char msg[1000];
    // printf("'%s', ", List[i]);
    SPRINTF(msg, "%.900s '%.50s',", dummy, List[i]);    
    STRCPY(dummy, msg);
  }
  //printf("'%s'\n ", List[i]);
  	  
  RFERROR2("%.900s and '%.50s'.", dummy, List[i]);  
 
 ErrorHandling:
  if (defaultvalue >= 0) {
    ans[0] = defaultvalue;
    for (k=1; k<maxlen_ans; k++) ans[k] = endvalue;
    return;
  }
  
  RFERROR1("'%.50s': no value given.", name);
}

int GetName(SEXP el, char *name, const char * List[], int n,
	    int defaultvalue) {
  int i;
  GetName(el, name, List, n, defaultvalue, defaultvalue, &i, 1);
  return i;
}


int GetName(SEXP el, char *name, const char * List[], int n) {
 return GetName(el, name, List, n, -1);
}


double ownround(double x) { return TRUNC((x + SIGN(x) * 0.5)); }


double lonmod(double x, double modulus) {  
  double 
    halfmodulus = 0.5 * modulus,
    y = x + modulus + halfmodulus;
  return Mod(y, modulus) - halfmodulus;
}
