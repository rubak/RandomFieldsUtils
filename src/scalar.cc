
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of system specific auxiliary functions

 Copyright (C) 2001 -- 2015 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/*

Makefile must be:

PKG_LIBS =  $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)  -march=native  -mssse3 
PKG_CXXFLAGS =  $(SHLIB_OPENMP_CXXFLAGS)  -march=native -mssse3 

 */



#define SIMD_AVAILABLE 1



#include "RandomFieldsUtils.h"
#include "General_utils.h"
#ifdef XXXSCHLATHERS_MACHINE

#ifdef SIMD_AVAILABLE
#include <x86intrin.h>
#endif
 


#include "kleinkram.h"

#define Nmodi 9
name_type modi = { "1x1", "2x2", "4x4", "8x8", "near", "simple", "precise", "kahan", "1x1p"};



typedef unsigned int uint32;


#define size 8
#define vectorlen (256 / (size * 8))
#define repet 8
#define VECTOR _mm256_loadu_pd
#define SET_0(NR) sum##NR = _mm256_setzero_pd()
#define P_0(NR) prod##NR = _mm256_setzero_pd()
#define SUMUP(NR, nr) sum##NR = _mm256_add_pd(sum##NR, sum##nr)
#define ADDF(NR) \
  sum##NR = _mm256_fmadd_pd(VECTOR(x + i + NR * vectorlen),\
			    VECTOR(y + i + NR * vectorlen), sum##NR)
#define ADDN(NR)							\
  prod##NR = _mm256_mul_pd(VECTOR(x + i + NR * vectorlen),		\
			   VECTOR(y + i + NR * vectorlen));		\
  sum##NR = _mm256_add_pd(sum##NR, prod##NR) 


#ifdef SIMD_AVAILABLE

#ifdef FMA_AVAILABLE
double avx_scalarproductDfma(double * x, double * y, int len) {
  int i = 0,
     lenM = len - (repet * vectorlen - 1);  
   __m256d SET_0(0);
   double *D  = (double *) &sum0;

  if ( len >= vectorlen * repet) {
    __m256d SET_0(1), SET_0(2), SET_0(3), SET_0(4), SET_0(5), SET_0(6),SET_0(7);
#if (7 != repet - 1)
  wrong repet length
#endif
   for (; i < lenM; i += repet * vectorlen) { 
     ADDF(0); ADDF(1); ADDF(2); ADDF(3); ADDF(4); ADDF(5); ADDF(6); ADDF(7); 
#if (7 != repet - 1)
  wrong repet length
#endif
    }
  SUMUP(0, 1); SUMUP(2, 3); SUMUP(4, 5); SUMUP(6, 7); SUMUP(0, 2); SUMUP(4, 6); SUMUP(0, 4);
#if (7 != repet - 1)
  wrong repet length
#endif
  }
  lenM = len - vectorlen + 1;
  for (; i < lenM; i += vectorlen) { // could unroll further
    ADDF(0);
  }
  double sum = D[0] + D[1] + D[2] + D[3];
#if (3 != vectorlen - 1)
  wrong vector length
#endif

  for (; i < len; ++i) sum += x[i] * y[i];
  return sum;
}
#endif


double avx_scalarproductDnearfma(double * x, double * y, int len) {
  // deutlich genauer zum 0 tarif
  int i = 0,
     lenM = len - (repet * vectorlen - 1);  
  __m256d SET_0(0), SET_0(1), SET_0(2), SET_0(3), SET_0(4), SET_0(5), SET_0(6),SET_0(7),
    P_0(0), P_0(1), P_0(2), P_0(3), P_0(4), P_0(5), P_0(6),P_0(7);

   double *D  = (double *) &sum0;

  if ( len >= vectorlen * repet) {
    for (; i < lenM; i += repet*vectorlen) {
      //
      ADDN(0); ADDN(1); ADDN(2); ADDN(3); ADDN(4); ADDN(5); ADDN(6); ADDN(7); 
 #if (7 != repet - 1)
  wrong repet length
#endif
    }
    SUMUP(0, 1); SUMUP(2, 3); SUMUP(4, 5); SUMUP(6, 7); SUMUP(0, 2); SUMUP(4, 6); SUMUP(0, 4);
    // SUMUP(0, 1); SUMUP(0, 2); SUMUP(0, 3); SUMUP(0, 4); SUMUP(0, 5); SUMUP(0, 6); SUMUP(0, 7);
#if (7 != repet - 1)
  wrong repet length
#endif
  }
  lenM = len - vectorlen + 1;
  for (; i < lenM; i += vectorlen) { // could unroll further
    ADDN(0);
  }

  double sum = D[0] + D[1] + D[2] + D[3];
#if (3 != vectorlen - 1)
  wrong vector length
#endif

  for (; i < len; ++i) sum += x[i] * y[i];
  return sum;
}
 

#define ADD(NR)								\
  prod0 = _mm256_mul_pd(VECTOR(x + i + NR * vectorlen),		\
			   VECTOR(y + i + NR * vectorlen));		\
  sum0 = _mm256_add_pd(sum0, prod0)
double avx_scalarproductD(double * x, double * y, int len) {
  int i = 0,
     lenM = len - (repet * vectorlen - 1);  
  __m256d SET_0(0), P_0(0);
   double *D  = (double *) &sum0;

  if ( len >= vectorlen * repet) {
  for (; i < lenM; i += repet*vectorlen) {
    //
    ADD(0); ADD(1); ADD(2); ADD(3); ADD(4); ADD(5); ADD(6); ADD(7); 
 #if (7 != repet - 1)
  wrong repet length
#endif
    }
  }

  lenM = len - vectorlen + 1;
  for (; i < lenM; i += vectorlen) { // could unroll further
    ADD(0);
  }
  
  double sum = D[0] + D[1] + D[2] + D[3];
#if (3 != vectorlen - 1)
  wrong vector length
#endif

  for (; i < len; ++i) sum += x[i] * y[i];
  return sum;
}


double avx_scalarproductDP(double * x, double * y, int len) {
  int i = 0,
     lenM = len - (repet * vectorlen - 1);  
  __m256d SET_0(0), SET_0(1), P_0(0);
   double *D  = (double *) &sum1;

  if ( len >= vectorlen * repet) {
    
    for (; i < lenM; ) {
      int lenMM = i + vectorlen * (repet * 10 + 1);
      if (lenMM > lenM) lenMM = lenM;
      sum0 = _mm256_mul_pd(VECTOR(x + i), VECTOR(y + i));
      i += vectorlen;
      for (; i < lenMM; i += repet*vectorlen) {
	ADD(0); ADD(1); ADD(2); ADD(3); ADD(4); ADD(5); ADD(6); ADD(7); 
#if (7 != repet - 1)
	wrong repet length
#endif
	  }
      sum1 = _mm256_add_pd(sum0, sum1);
    }
  }
  
 lenM = len - vectorlen + 1;
 for (; i < lenM; i += vectorlen) { // could unroll further
    prod0 = _mm256_mul_pd(VECTOR(x + i), VECTOR(y + i));
    sum1 = _mm256_add_pd(sum1, prod0);
  }
  
  double sum = D[0] + D[1] + D[2] + D[3];
#if (3 != vectorlen - 1)
  wrong vector length
#endif

    for (; i < len; ++i) {
      // printf("final %d\n", i);
      sum += x[i] * y[i];
    }
  return sum;
}




#define ADDK(NR)								\
  prod0 = _mm256_mul_pd(VECTOR(x + i + NR * vectorlen),		\
			   VECTOR(y + i + NR * vectorlen));		\
  sum2 = _mm256_sub_pd(prod0, sum1);\
  sum3 = _mm256_add_pd(sum0, sum2);		\
  sum1 = _mm256_sub_pd(sum3, sum0);		\
  sum0 = sum3;					\
  sum1 = _mm256_sub_pd(sum1, sum2);
 double avx_scalarproductDK(double * x, double * y, int len) {
   // Kahan enhanced
  int i = 0,
     lenM = len - (repet * vectorlen - 1);  
  __m256d SET_0(0), // sum
    SET_0(1),  // c
    SET_0(2), // y
    SET_0(3),  // t
    P_0(0),
    P_0(1);
   double *D  = (double *) &sum0;

  if ( len >= vectorlen * repet) {
  for (; i < lenM; i += repet*vectorlen) {
    //
    ADDK(0); ADDK(1); ADDK(2); ADDK(3); ADDK(4); ADDK(5); ADDK(6); ADDK(7);
 #if (7 != repet - 1)
  wrong repet length
#endif
    }
  }
 lenM = len - vectorlen + 1;
 for (; i < lenM; i += vectorlen) { // could unroll further
    ADDK(0);
 }
 sum0 = _mm256_add_pd(sum0, prod1); 
  
  double sum = D[0] + D[1] + D[2] + D[3];
#if (3 != vectorlen - 1)
  wrong vector length
#endif

  for (; i < len; ++i) sum += x[i] * y[i];
  return sum;
}

// end if simd
#endif
 
 
double scalarproductf64( double * v1,  double * v2, int N){
  double *endv1 = v1 + N,
    sum = 0;
    for(; v1!= endv1; v1++, v2++) sum+= v2[0] * v1[0];
    return sum;
}
 
 
 
double scalarproductf64P( double * v1,  double * v2, int N){
  double //*endv1 = v1 + N,
    sum = 0;
#ifdef DO_PARALLEL
#pragma omp parallel for reduction(+:sum)
#else
  ERR("parallel not allowed");
#endif
  for(int i=0; i<=N; i++) sum += v2[i] * v1[i];
    return sum;
}
 
 
double scalarproduct2by2f64( double * v1,  double * v2, int N){
    double *endv1 = v1 + N,
      sum = 0;
    for(; v1!= endv1; v1+=2, v2+=2) {
        sum+= v2[0] * v1[0] + v2[1] * v1[1];
    }
    return sum;
}
 
 
double scalarproduct4by4f64( double * v1,  double * v2, int N){
    double*endv1 = v1 + N,
      sum = 0;
    for(; v1 < endv1; v1+=4, v2+=4) {
      sum+= v2[0] * v1[0] + v2[1] * v1[1] + v2[2] * v1[2]+ v2[3] * v1[3];
    }
    return sum;
}

 
double scalarproduct8by8f64( double * v1,  double * v2, int N){
    double *endv1 = v1 + N,
      sum = 0;
    for(; v1!= endv1; v1+=8, v2+=8) {
      sum+= v2[0] * v1[0] + v2[1] * v1[1]+ v2[2] * v1[2] + v2[3] * v1[3] +
	v2[4] * v1[4] + v2[5] * v1[5]+ v2[6] * v1[6]+ v2[7] * v1[7];
    }
    return sum;
}
 

SEXP scalarX(SEXP x, SEXP y, SEXP mode) {
  int len = length(x);
  if (length(y) != len) ERR("x and y differ in length");
  int n = Match((char*) CHAR(STRING_ELT(mode, 0)), modi, Nmodi);
  if (n < 0) ERR("unknown modus");
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, 1));
  double *ans = REAL(Ans);
  switch(n) {
  case 0 : *ans = scalarproductf64(REAL(x), REAL(y), len); break;
  case 1 : *ans = scalarproduct2by2f64(REAL(x), REAL(y), len); break;
  case 2 : *ans = scalarproduct4by4f64(REAL(x), REAL(y), len); break;
  case 3 : *ans = scalarproduct8by8f64(REAL(x), REAL(y), len); break;
  case 4 :
#ifdef SIMD_AVAILABLE
    *ans = avx_scalarproductDnearfma(REAL(x), REAL(y), len); break;
    #else
    BUG;
    #endif
  case 5 :
    #ifdef SIMD_AVAILABLE
    *ans = avx_scalarproductD(REAL(x), REAL(y), len); break;
     #else
    BUG;
   #endif
  case 6 :
    #ifdef SIMD_AVAILABLE
    *ans = avx_scalarproductDP(REAL(x), REAL(y), len); break;
    #else
    BUG;
    #endif
  case 7 :
    #ifdef SIMD_AVAILABLE
    *ans = avx_scalarproductDK(REAL(x), REAL(y), len); break;
    #else
    BUG;
    #endif
  case 8 : *ans = scalarproductf64P(REAL(x), REAL(y), len); break;
  default : BUG;
  }
 
  UNPROTECT(1);
  return Ans;
}
  
#else
SEXP scalarX(SEXP x, SEXP y, SEXP mode) { BUG; }


#endif // SCHLATHERS_MACHINE
