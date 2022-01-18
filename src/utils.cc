
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of system specific auxiliary functions

 Copyright (C) 2001 -- 2021 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.
g
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include "Basic_utils_local.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include "RandomFieldsUtils.h"
#include "kleinkram.h"
#include "zzz_RandomFieldsUtils.h"
#include "Utils.h"
#include "xport_import.h"
#include "extern.h"


AVAILABLE_SIMD


double *ToRealI(SEXP X, bool *create) {
  KEY_type *KT = KEYT();
  if (TYPEOF(X) == REALSXP) { 
    *create = false;
    return REAL(X);
  }
  // TO DO !!
  //HELPINFO("Better use 'double' as storage mode (for one of the arguments).");
  int len = length(X); 
  double *y;
  if (create || KT->ToRealN < len) {
    y = (double *) MALLOC(sizeof(double) * len);
    if (y == NULL) ERR1("not enough memory for an %d vector of doubles", len);
    if (!create) {
      FREE(KT->ToRealDummy);
      KT->ToRealDummy = y;
      KT->ToRealN = len;
    }
  } else y = KT->ToRealDummy;
  int *x;
  if (TYPEOF(X)==INTSXP) x=INTEGER(X); else x=LOGICAL(X);
  for (int i=0; i<len; i++) y[i] = (double) x[i];
  return y;
}

double *ToReal(SEXP X) {
  bool ToFalse[1] = { false };
 if (TYPEOF(X) == REALSXP) return REAL(X);
  return ToRealI(X, ToFalse);
}


#ifdef __cplusplus
// OBSOLETE_RFU
extern "C" {
#endif

int *ToIntI(SEXP X, bool *create, bool round) {
   KEY_type *KT = KEYT();
  if (TYPEOF(X) == INTSXP) {
    *create = false;
    return INTEGER(X);
  }
  if (TYPEOF(X) == LGLSXP) {
    *create = false;
    return LOGICAL(X);
  }
  int len = length(X);
  // TO DO !!
  //  if (len > 100 || PL > 1)
  //    HELPINFO("Better use 'integer' as storage mode (for one of the arguments).");
  int *y;
  if (*create || KT->ToIntN < len) {
    y = (int *) MALLOC(sizeof(int) * len);    
    if (y == NULL) ERR1("not enough memory for an %d vector of integers", len);
    if (!*create) {
      FREE(KT->ToIntDummy);
      KT->ToIntDummy = y;
      KT->ToIntN = len;
    }
  } else y = KT->ToIntDummy;
  double *x = (double *) REAL(X);
  if (round) for (int i=0; i<len; i++) y[i] = (int) ROUND(x[i]);
  else for (int i=0; i<len; i++) y[i] = (int) x[i];
  return y;
}
#ifdef __cplusplus
}
#endif
  

int *ToInt(SEXP X) {
  bool ToFalse[1] = { false };
  return ToIntI(X, ToFalse, false);
}




SEXP DivByRow(SEXP M, SEXP V) {
  Long
    l = length(V),
    r = nrows(M),
    c = ncols(M);

  double *m = REAL(M),
    *v = REAL(V);
  
  if (l != c) ERR0("vector does not match matrix");
  for (Long j=0; j<c; j++) {
    double vj = v[j];
    for (Long i=0; i<r; i++) {
      *(m++) /= vj;
    }
  }

  return M;
}

#define algn_general(X)  ((1U + (uintptr_t) (((uintptr_t) X - 1U) / BytesPerBlock)) * BytesPerBlock)

double static inline *algn(double *X) {
  assert(algn_general(X)>=(uintptr_t)X); return (double *) algn_general(X);
}

#if defined SSE41 || defined AVX2
int static inline *algnInt(int *X) {
  assert(algn_general(X)>=(uintptr_t)X); return (int *) algn_general(X);
}
#endif



void colMaxsIint(int *M, Long r, Long c, int *ans) {
  if (r < 32
#if defined AVX2
      || !avx2Avail
#elif defined  SSE41
      || !sse41Avail
#endif      
       ) {
    for (Long i=0; i<c; i++) {
      int *m = M + r * i,
	dummy = m[0];    
      for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
      ans[i] = dummy;
    }
    return;
  }
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(static)
#endif    
  for (Long i=0; i<c; i++) {
     int dummy,
      *m = M + r * i;
#if defined SSE41 || defined AVX2    
     int *start = algnInt(m),
       *end = m + r;
    uintptr_t End = (uintptr_t) (end - integers);
    if ((uintptr_t) start < End) {
      BlockType *m0 = (BlockType0*) start,
	Dummy = LOAD((BlockType0*) m0);
      for (m0++ ; (uintptr_t) m0 < End; m0++) {
	Dummy = MAXINTEGER(Dummy, LOAD(m0));
      }
      int *d = (int *) &Dummy;
      dummy = d[0];
      dummy = MAX(dummy, d[1]);
      dummy = MAX(dummy, d[2]);
      dummy = MAX(dummy, d[3]);
#if defined AVX2
      dummy = MAX(dummy, d[4]);
      dummy = MAX(dummy, d[5]);
      dummy = MAX(dummy, d[6]);
      dummy = MAX(dummy, d[7]);
#endif // AVX
      for ( ; m<start; m++) dummy = MAX(dummy, *m);
      m = (int *) m0;
      for ( ; m<end; m++) dummy = MAX(dummy, *m);
    } else {
      dummy = m[0];    
      for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
    }
#else // not SSE4
    dummy = m[0];    
    for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
#endif    
    ans[i] = dummy;
  }
}


void colMaxsI(double *M, Long r, Long c, double *ans) {
  if (r < 16
#if defined AVX2
      || !avx2Avail
#elif defined  SSE2
      || !sse2Avail
#endif      
      ) {
    for (Long i=0; i<c; i++) {
      double *m = M + r * i,
	dummy = m[0];    
      for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
      ans[i] = dummy;
    }
    return;
  }
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(static)
#endif  
  for (Long i=0; i<c; i++) {
    double dummy,
      *m = M + r * i;
#if defined SSE2
    double *start = algn(m),
      *end = m + r;
    uintptr_t End = (uintptr_t) (end - doubles);
    if ((uintptr_t) start < End) {
      Double * m0 = (Double*) start,
	Dummy = (Double) LOAD((BlockType0*) m0);
      for (m0++ ; (uintptr_t) m0 < End; m0++) {
	Dummy = MAXDOUBLE(Dummy, (Double) LOAD((BlockType0*) m0));
      }
      double *d = (double *) &Dummy;
      dummy = d[0];
      dummy = MAX(dummy, d[1]);
#if defined AVX
      dummy = MAX(dummy, d[2]);
      dummy = MAX(dummy, d[3]);
#endif
      for ( ; m<start; m++) dummy = MAX(dummy, *m);
      m = (double *) m0;
      for ( ; m<end; m++) dummy = MAX(dummy, *m);
    } else {
      dummy = m[0];    
      for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
    }
#else
    dummy = m[0];    
    for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
#endif    
    ans[i] = dummy;
  }
}


void colMaxsI256(double *M, Long r, Long c, double *ans);
void colMaxsIint256(int *M, Long r, Long c, int *ans);
SEXP colMaxs(SEXP M) {
  Long
    r = nrows(M),
    c = ncols(M);
  if (r == 0) return R_NilValue;
  SEXP Ans;
  if (TYPEOF(M) == REALSXP) {
    PROTECT(Ans = allocVector(REALSXP, c));
    if (avx2Avail) colMaxsI256(REAL(M), r, c, REAL(Ans));
     else colMaxsI(REAL(M), r, c, REAL(Ans));
  } else {
    bool i = TYPEOF(M) == INTSXP;
    PROTECT(Ans = allocVector(i ? INTSXP : LGLSXP, c));
    int *m, *a;
    if (i) {
      m = INTEGER(M);
      a = INTEGER(Ans);
    } else {
      m = LOGICAL(M);
      a = LOGICAL(Ans);
    }
    if (avx2Avail) colMaxsIint(m, r, c, a);
    else colMaxsIint256(m, r, c, a);
  }
  UNPROTECT(1);
  return Ans;
}


SEXP rowProd(SEXP M) {
  Long
    r = nrows(M),
    r4 = r / 4,
    c = ncols(M);
  if (r == 0) return R_NilValue;
  SEXP Ans;
  if (TYPEOF(M) == REALSXP) {
    PROTECT(Ans = allocVector(REALSXP, r));
    double *ans = REAL(Ans),
      *m = REAL(M);
    MEMCOPY(ans, m, sizeof(double) * r);
    m += r;
    for (Long ic=1; ic<c; ic++) {
      double *a = ans;
      for (Long ir=0; ir<r4; ir++) {
	*(a++) *= *(m++);
	*(a++) *= *(m++);
	*(a++) *= *(m++);
	*(a++) *= *(m++);
      }
      for (Long ir=r4 * 4; ir<r; ir++) *(a++) *= *(m++);
    }
  } else {
    // printf("type = %d", TYPEOF(M));
    RFERROR("transform to double first") ;
  }
  UNPROTECT(1);
  return Ans;
}

SEXP rowMeansX(SEXP M, SEXP Weight) {
  // todo : SSE2 / AVX
  Long
    r = nrows(M),
    c = ncols(M);
  if (r == 0 || c == 0) return R_NilValue;
  if (length(Weight) != c && length(Weight) != 0)
    ERR0("Length of 'weight' must equal number of columns of 'x'.");
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, r));
  double *ans = REAL(Ans);
  for (Long j=0; j<r; j++) ans[j] = 0.0;
  if (length(Weight) == 0) {    
#define for1					\
    for (Long i=0; i<c; i++, m+=r) {			\
      for (Long j=0; j<r; j++) ans[j] += (double) m[j];	\
    }
  
    if (TYPEOF(M) == REALSXP) { double *m = REAL(M); for1; }
    else {
      int *m;
      if (TYPEOF(M) == INTSXP) m = INTEGER(M); else m = LOGICAL(M);
      for1;
    }
    
  } else {    
    double *weight = ToReal(Weight);
#define for2							\
    for (Long i=0; i<c; i++, m+=r) {				\
      double dummy = weight[i]; /* load1(weight); MULTDOUBLE */ \
      for (Long j=0; j<r; j++) ans[j] += (double) m[j] * dummy;	\
    }

    if (TYPEOF(M) == REALSXP) { double *m = REAL(M); for2; }
    else {
      int *m;
      if (TYPEOF(M) == INTSXP) m = INTEGER(M); else m = LOGICAL(M);
      for2;
    }
    
    if (TYPEOF(Weight) != REALSXP) { FREE(weight); }
  }
  double invc = 1.0 / (double) c;
  for (Long j=0; j<r; j++) ans[j] *= invc;
  UNPROTECT(1);
  return Ans;
}

SEXP dbinorm(SEXP X, SEXP Sigma) { // 12'41
  Long nrow,
    ncol = 2;
  double *x, *y;
  if (TYPEOF(X) == VECSXP) {
    if (length(X) != ncol) BUG;
    SEXP xx = VECTOR_ELT(X, 0);
    nrow = length(xx);
    x = REAL(xx);
    y = REAL(VECTOR_ELT(X, 1));
  } else {
    if (isMatrix(X)) {
      if (ncols(X) != ncol) BUG;
      nrow = nrows(X);
    } else if (isVector(X)) {
      if (length(X) != ncol) BUG;
      nrow = 1;
    } else BUG;
    x = REAL(X);
    y = x + nrow;
  }

  
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, nrow));
  double *ans = REAL(Ans);
  //  Long nrow4 = nrow - 4;
  if (length(Sigma) == 0) {
    double invtwopi = 1.0 / TWOPI;
    /*
      minushalfX[4] ={-0.5, -0.5, -0.5, -0.5},
      invtwopiX [4] = {invtwopi, invtwopi, invtwopi, invtwopi};
    Long i=0;

#define atonce 4
    __m256d minushalf4 = LOADuDOUBLE(minushalfX),
       invtwopi4 = LOADuDOUBLE(invtwopiX);
      
    for (; i<nrow4; i+=atonce) {
      __m256d x4 = LOADuDOUBLE(x + i);
      double *xx4 = (double *) &x4;
      x4 = MULTDOUBLE(x4, x4);
      {
	__m256d y4 = LOADuDOUBLE(y + i);
	y4 = MULTDOUBLE(y4, y4);
	x4 = ADDDOUBLE(x4, y4);
      }
      x4 = MULTDOUBLE(minushalf4, x4);
      xx4[0] = EXP(xx4[0]);
      xx4[1] = EXP(xx4[1]);
      xx4[2] = EXP(xx4[2]);
      xx4[3] = EXP(xx4[3]);
      x4 = MULTDOUBLE(x4, invtwopi4);
      STOREuDOUBLE(ans + i, x4);
    }
    */
    for (Long i=0; i<nrow; i++) 
      ans[i] = EXP(-0.5 * (x[i] * x[i] + y[i] * y[i])) * invtwopi;
    } else {
    double *sigma=REAL(Sigma),
      sigma1 = sigma[0],
      sigma4 = sigma[3],
      inv2piSrtS = 1.0 / (TWOPI * SQRT(sigma1 * sigma4)),
      invS1half = 0.5 / sigma1,
      invS4half = 0.5 / sigma4;
    if (sigma[1] == 0.0 && sigma[2] == 0.0) {
      for (Long i=0 ; i<nrow; i++)
	ans[i] = EXP(- (x[i] * x[i] * invS1half + y[i] * y[i] * invS4half) )
	  * inv2piSrtS;
    } else BUG;
  }
  UNPROTECT(1);
  return Ans;
}




SEXP test(SEXP AA, SEXP CC, SEXP X) {
  KEY_type *KT = KEYT();
  int cores = KT->global_utils.basic.cores;
  Long nrow = nrows(AA),
    ncol = ncols(AA),
    dim = length(X),
    k = MIN(ncol / 2, nrow),
    m = MAX(ncol, nrow);
  
  double
    eps = 1e-14,
    *A = REAL(AA),
    *C = REAL(CC),
    *x = REAL(X),    
    z[2],
    *a[2] = {(double*) MALLOC(m * m * sizeof(double)),
	       (double*) MALLOC(m * m * sizeof(double))};

  if (ncols(CC) != nrows(CC) ||  ncols(CC) != ncol) BUG;
  if (length(X) != nrow) BUG;

  for (int i=0; i<=17; i++) {
    for (int j=0; j<=1; j++) {
      extern bool obsolete_package_in_use;
      SetLaMode(j == 0 || obsolete_package_in_use
		? LA_INTERN : LA_R, cores);
      switch(i) {
      case 1: z[j] = XkCXtl(A, C, nrow, ncol, nrow / 3, nrow / 4, cores); break;
      case 2: XCXt(A, C, a[j], nrow, ncol, cores); break;
      case 3: AtA(A, nrow, ncol, a[j], cores); break;
      case 4: xA(x, A, nrow, ncol, a[j], cores); break;
      case 5: xA_noomp(x, A, nrow, ncol, a[j]); break;
	//    case : xA(x1, x2,  A, nrow, ncol, a[j]1,  a[j]2); break;
      case 6: z[j] = xAx(x, C, nrow, cores); break;
      case 7: Ax(A, C, nrow, ncol, a[j], cores); break;// C genuegend lang. Reicht.
      //    case 8: Ax(A, x, x2, nrow, ncol, a[j]1,  a[j]2); break;
      case 8: z[j] =xUy(x, C, A, dim, cores); break; // A genuegend lang. Reicht.
      case 9: z[j] =xUxz(x, C, dim, a[j], cores); break;
      case 10: z[j] =x_UxPz(x, C, A, dim,cores); break; // A genuegend lang. Reicht.
      case 11: z[j] =xUx(x, C, dim, cores); break;
      case 12: matmult(A, C, a[j], nrow, ncol, k, cores); break;
      case 13: matmulttransposed(A, C, a[j], ncol, nrow, k, cores); break;
	//case : matmulttransposedInt(int *A, int *B, int *c, ncol, ncol, k); break; 
      case 14: matmult_2ndtransp(A, C, a[j], nrow, ncol, k, cores); break;
      case 15: matmult_2ndtransp(A, C, a[j], nrow, ncol, cores); break;
      case 16: matmult_tt(A, C, a[j], ncol, nrow, k,cores); break;
	//     case 17: z[j]=  scalar(A, C, ncol); break;	
      default: BUG;
      }

      int size = 0;
      switch(i) {
      case 1: case 6: case 8: case 9:case 10: case 11: case 17:
	if (FABS(z[0] - z[1])> eps) { PRINTF("i=%d", i); BUG; }
	break;
      case 2:  size = ncol * ncol;
	break;
      case 3: case 15: size = nrow * nrow;
	break;
      case 4 : case 5: size = ncol;
	break;
      case 7 : size = nrow;
	break;
      case 12: case 13: case 14: case 16: size = nrow * k;
	break;
      default: BUG;
      }
      for (int p=0; p<size; p++)
	if (FABS(a[0][p] - a[1][p]) > eps)  { PRINTF("i=%d, %d", i, p); BUG; }
    }
  }

 FREE(a[0]);
 FREE(a[1]);
  
  return R_NilValue;
}



SEXP quadratic(SEXP A, SEXP x) {
  KEY_type *KT = KEYT();
  int cores = KT->global_utils.basic.cores;
  SEXP ans;
  int len = length(x);
  if (len != nrows(A) || len != ncols(A)) ERR0("'x' and 'A' do not match.");
  PROTECT(ans = allocVector(REALSXP, 1));
  REAL(ans)[0] = xAx(REAL(x), REAL(A), len, cores);
  UNPROTECT(1);
  return ans;
}

SEXP dotXV(SEXP M, SEXP V) {
  Long
    r = nrows(M),
    c = ncols(M),
    l = length(V)
    ;
  if (l != r) ERR0("X and v do not match");
  if (r == 0) return R_NilValue;
  SEXP Ans;
  PROTECT(Ans = allocMatrix(REALSXP, r, c));

  // bringt nix
  //#ifdef DO_PARALLEL
  //#p ragma omp parallel for num_threads(CORES) 
  //#endif  
  for (Long i=0; i<c; i++) {
    //  printf("i=%d\n", i);
#if defined SSE2_DONOTUSE_AS_SLOWER
    double 
      *ans = REAL(Ans) + r * i,
      *v = REAL(V),
      *m = REAL(M) + r * i,
      *end = m + r - doubles;
    for ( ; m < end; m += doubles, ans += doubles, v += doubles)
      STOREuDOUBLE(ans, MULTDOUBLE(LOADuDOUBLE(m), LOADuDOUBLE(v)));
    end += doubles;
    for (; m < end; m++) *ans = *m * *v;
#else
    double
      *ans = REAL(Ans) + r * i,
      *v = REAL(V),
      *m = REAL(M) + r * i;
    for (Long j=0; j<r; j++) {
      ans[j] = m[j] * v[j];
    }
     
#endif    
  }

  UNPROTECT(1);
  return Ans;
}




SEXP debuggingLevel() {
  SEXP ans;
  PROTECT(ans = allocVector(INTSXP, 1));
#ifdef SCHLATHERS_MACHINE
  INTEGER(ans)[0] = 1;
#else
  INTEGER(ans)[0] = 0;
#endif  
  UNPROTECT(1);
  return ans;
}
// for debugging only
SEXP DebugCall() {
  //  return R_NilValue;
  //   KEY_type *KT = KEYT();						
  //  assert((KT->n_data_names == 0) xor (KT->data_names != NULL)); 
  //  assert((KT->n_coord_names == 0) xor (KT->coord_names != NULL));
  //  assert((KT->n_data_idx == 0) xor (KT->data_idx != NULL));	
  //  assert((KT->n_coord_idx == 0) xor (KT->coord_idx != NULL));
  return R_NilValue;
}






#define Nmodi 9
name_type modi = { "1x1", "2x2", "4x4", "8x8", "near", "simple", "precise", "kahan", "1x1p"};

 
double scalarprod( double * v1, double * v2, Long N){
  double *endv1 = v1 + N,
    sum = 0;
  for(; v1!= endv1; v1++, v2++) sum +=  v2[0] * v1[0];
  return sum;
}
 
 
double scalarprod2by2( double * v1, double * v2, Long N){
  double *endv1 = v1 + (N / 2) * 2,
    *end = v1 + N,
    sum = 0;
  for(; v1 < endv1; v1 += 2, v2 += 2) sum += v2[0] * v1[0] + v2[1] * v1[1];
  if (v1 < end) sum += v2[0] * v1[0]; 
  return sum;
}
 
 
double scalarprod4by4( double * v1, double * v2, Long N){
  // printf("4by4 %d %d %d\n", sse, sse2, avx);
  double*endv1 = v1 + (N / 4) * 4,
    *end = v1 + N,
    sum = 0;
  for(; v1 < endv1; v1 += 4, v2 += 4)
    sum += v2[0] * v1[0] + v2[1] * v1[1] + v2[2] * v1[2]+ v2[3] * v1[3];
  for(; v1 < end; v1++, v2++) sum += v2[0] * v1[0];        
  return sum;
}

 
double scalarprod8by8( double * v1, double * v2, Long N){
  double
    *endv1 = v1 + (N / 8) * 8,
    *end = v1 + N,
    sum = 0.0;
  for(; v1 < endv1; v1 += 8, v2 += 8)
    sum += v2[0] * v1[0] + v2[1] * v1[1]+ v2[2] * v1[2] + v2[3] * v1[3] +
      v2[4] * v1[4] + v2[5] * v1[5]+ v2[6] * v1[6]+ v2[7] * v1[7];
  for(; v1 < end; v1++, v2++) sum +=  v2[0] * v1[0];        
  return sum;
}


void avx_scalarprodM(double * x, double * y, Long len, double *res);
double avx_scalarprodDnearfma(double * x, double * y, Long len);
double avx_scalarprodD(double * x, double * y, Long L);
double avx_scalarprodDopt(double * x, double * y, Long L);
double avx_scalarprodDP(double * x, double * y, Long L) ;
double avx_scalarprodDK(double * x, double * y, Long L);
double scalarX(double *x, double *y, Long len, Long n) {
  // parallel lohnt i.A. nicht: 28.2.20121 alles was parallel ist, rausgeworfen
  assert(n >= 0);
  //  __m128 a, b;  a = _mm_add_ps ((__m128) a, (__m128) b);
  //   __m128i c, d; c = _mm_add_epi16 ((__m128i) c, (__m128i) d);
  //  __m256d e, f;   e = _mm256_add_pd ( e,  f);
   //  printf("n=%d %d ",n, avx);
  switch(n) {
    //  printf("%d\n", n);
  case SCALAR_AVX :
    //  printf(" # %d ", avx);
    if (avxAvail) return avx_scalarprodD(x, y, len); // best one kernel
    break;
  case 2 : return scalarprod(x, y, len);
  case 3 : return scalarprod2by2(x, y, len); 
  case 4 : return scalarprod8by8(x, y, len); 
    //  case 5 :
    //#ifdef FMA_AVAILABLE
    //   return avx_scalarprodDfma(x, y, len);
    //#endif    
  case SCALAR_NEARFMA :
    if (avxAvail) return avx_scalarprodDnearfma(x, y, len);
    break;
  case 7 :
    if (avxAvail) return avx_scalarprodDP(x, y, len);  //best
    break;
  case SCALAR_KAHAN :
    if (avxAvail) return avx_scalarprodDK(x, y, len); // kahan   
    break;

    /*
  case 10:
    if (avx) return avx_scalarprodDopt(x, y, len); // best one kernel
    break;

  case 11:
    double result[6];
    if (avx) {
      avx_scalarprodM(x, y, len, result); // best one kernel      
      return result[0];
    }
    break;
    */

  case SCALAR_BASE :
  default :
    {}
  }
  return scalarprod4by4(x, y, len);
}
  


SEXP scalarR(SEXP x, SEXP y, SEXP Mode) { // unused
  Long len = length(x);
  if (length(y) != len) ERR0("x and y differ in length");
  int mode;
  if (length(Mode) == 0) mode = -1;
  else if (INTSXP==TYPEOF(Mode)) mode = INTEGER(Mode)[0];
  else mode = Match((char*) CHAR(STRING_ELT(Mode, 0)), modi, Nmodi);
  SEXP Ans;

  if (isMatrix(x)) {
    Long nc = ncols(x);
    PROTECT(Ans = allocVector(REALSXP, nc * (nc - 1) / 2));
    double *ans = REAL(Ans);
    *ans = scalarX(REAL(x), REAL(y), len, 11); // no PROTECT( needed
    UNPROTECT(1);
  } else {
    PROTECT(Ans = allocVector(REALSXP, 1));
    double *ans = REAL(Ans);
    *ans = scalarX(REAL(x), REAL(y), len, mode); // no PROTECT( needed
    UNPROTECT(1);
  }
  return Ans;
}


SEXP crossprodX(SEXP X, SEXP Y, SEXP mode) {
  KEY_type *KT = KEYT();
  int cores = KT->global_utils.basic.cores;
  Long n, nrow,
    len,
    lenY,
    ncol;
  if (isMatrix(X)) {
    nrow = ncols(X);
    len = nrows(X);
  } else {
    nrow = 1;
    len = length(X);
  }
  if (isMatrix(Y)) {
    ncol = ncols(Y);
    lenY = nrows(Y);
  } else {
    ncol = 1;
    lenY = length(Y);
  }
  if (lenY != len) ERR0("sizes of 'x' and 'y' do not match");
  if (length(mode) == 0) n = SCALAR_DEFAULT;
  else {
    n = INTEGER(mode)[0];
    if (n < 0) n =  SCALAR_DEFAULT;
  }
  SEXP Ans; 
  PROTECT(Ans = allocMatrix(REALSXP, nrow, ncol));
  double *ans = REAL(Ans),
    *x = REAL(X),
    *y = REAL(Y);

  if (x == y) AtA(x, len, ncol, ans, cores);
  else matmulttransposed(x, y, ans, len, nrow, ncol, cores);

  UNPROTECT(1);
  return Ans;
}


void avx_linearprodD( double * v1,  double v2, Long N, double *inout);
void linearprod2by2( double * v1,  double v2, Long N, double *inout){
  double *endv1 = v1 + (N / 2) * 2,
    *end = v1 + N;
  for(; v1 < endv1; v1+=2, inout+=2) {
      inout[0] += v2 * v1[0];
      inout[1] += v2 * v1[1];
  }
  if (v1 < end) inout[0] += v2 * v1[0];
}
 

void linearX(double *x, double y, Long len, double *inout, Long n) {
  switch(n) {
  case LINEAR_AVX :
    if (avxAvail) { avx_linearprodD(x, y, len, inout); return; }
    break; 
  case LINEAR_BASE:
  default :
    {}
  }
  linearprod2by2(x, y, len, inout); 
}
  
