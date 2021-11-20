


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

#ifndef rfutils_init_H
#define rfutils_init_H 1


 /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
  */

#include "Utils.h"
#include "options.h"


#ifdef HAVE_VISIBILITY_ATTRIBUTE
  # define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
  # define attribute_hidden
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MY_PACKAGE "RandomFieldsUtils"
  //#define MY_ACRONYM XX
#include "zzz_calls.h"

  /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
!!!!! auch kein mit MALLOC kreiertes Objekt  !!!!
  */

#define AttachMessageN 2000
#define LEN_OPTIONNAME 201 // zwingend ungerade
typedef void (*setoptions_fctn) (int, int, SEXP, char[LEN_OPTIONNAME],
				 bool, bool);
typedef void (*getoptions_fctn) (SEXP, int, bool);
typedef void (*finalsetoptions_fctn) ();
typedef void (*deleteoptions_fctn) (bool);

  
#define MATERN_NU_THRES 100
#define BESSEL_NU_THRES 100
#define LOW_MATERN 1e-20
#define LOW_BESSEL 1e-20

  DECLARE1(void, utilsoption_DELETE, utilsoption_type *, S)
  DECLARE1(void, utilsoption_NULL, utilsoption_type *, S)
    
  DECLARE1(void, solve_DELETE, solve_storage**, S)
  DECLARE1(void, solve_NULL, solve_storage*, x)
  DECLARE7(int, solvePosDef, double*, M, int, size, bool, posdef, 
	   double *, rhs, Long, rhs_cols, double *, logdet, solve_storage *, PT)
  DECLARE8(int, solvePosDefSp, double *, M, int, size, bool, posdef,
	   double *, rhs, Long, rhs_cols, double *,logdet,
	   solve_storage *, Pt, solve_options *,sp)
  //  DECLARE8(int, solvePosDefResult, double*, M, int, size, bool, posdef, 
  //	   double *, rhs, int, rhs_cols, double *, result, double*, logdet, 
  //	   solve_storage*, PT)
  DECLARE4(int, sqrtPosDefFree, double *, M, int, size, solve_storage *, pt,
	   solve_options *, sp)
  DECLARE3(int, sqrtRHS, solve_storage *, pt, double*, RHS, double *, res)
  DECLARE2(int, invertMatrix, double *, M, int, size)
  DECLARE2(double, StruveH, double, x, double, nu)
  DECLARE3(double, StruveL, double, x, double, nu, bool, expScale1d)
  DECLARE1(double, I0mL0, double, x)
  DECLARE3(double, WM, double, x, double, nu, double, factor)
  DECLARE3(double, DWM, double, x, double, nu, double, factor)
  DECLARE3(double, DDWM, double, x, double, nu, double, factor)
  DECLARE3(double, D3WM, double, x, double, nu, double, factor)
  DECLARE3(double, D4WM, double, x, double, nu, double, factor)
  DECLARE4(double, logWM, double, x, double, nu1, double, nu2, double, factor)
  DECLARE1(double, Gauss, double, x)
  DECLARE1(double, DGauss, double, x)
  DECLARE1(double, DDGauss, double, x)
  DECLARE1(double, D3Gauss, double, x)
  DECLARE1(double, D4Gauss, double, x)
  DECLARE1(double, logGauss, double, x)
  DECLARE0(int, cores1)
  DECLARE0(int, cpus)
  
  DECLARE2(void, linkRFUoptions, utilsoption_type **, up, const char ***, R_T)
  DECLARE2(void, detachRFUoptions, const char **, prefixlist, int, N)
  DECLARE14(void, attachRFUoptions, char *, name,
	    const char **, prefixlist, int, N, 
	    const char ***, all, int *, allN, setoptions_fctn, set, 
	    finalsetoptions_fctn, final, getoptions_fctn, get,
	    deleteoptions_fctn, del, int, pl_offset, bool, basicopt,
	    install_modes, gpu_needs, Uint, avx_info, int, version)

  DECLARE3(void, sorting, double*, data, int, len, usr_bool, NAlast)
  DECLARE3(void, sortingInt, int*, data, int, len, usr_bool, NAlast)
  DECLARE4(void, ordering, double*, data, int, len, int, dim, int *, pos)
  DECLARE4(void, orderingInt, int*, data, int, len, int, dim, int *, pos)

    DECLARE4(double, scalarX, double *, x, double *, y, Long, len, Long, n)
  //  DECLARE4(int, scalarInt, int *, x, int *, y, int, len, int, n)
  DECLARE2(double, detPosDef, double *, M,  int, size) // destroys M!
  DECLARE3(double, detPosDefsp, double *, M, int, size, solve_options *, sp)
  DECLARE8(int, XCinvXdet,double*, M, int, size, double *,X, Long, X_cols,
	  double *, XCinvX, double *, det, bool, log, solve_storage, *PT)
  DECLARE10(int, XCinvYdet,double*, M, int, size, bool, posdef,
	    double *, X, double *, Y, Long, cols,
	    double *, XCinvY, double *, det, bool, log, solve_storage, *PT)
  //  DECLARE5(double, XCinvXlogdet, double *, M, int, size, double *, X,
  //	   int, X_cols, solve_storage *, PT)
  DECLARE2(bool, is_positive_definite, double *, C, Long, dim)
  DECLARE2(void, chol2inv, double *, MPT, int, size)
  DECLARE2(int, chol, double *, MPT, int, size)
  DECLARE1(void, pid, int *, i)
  DECLARE0(bool, parallel)
  DECLARE1(void, sleepMicro, int *, i)
  // DECLARE7(int, cholGPU, bool, copy, double*, M, int, size, double*, rhs, int, rhs_cols, double *, LogDet, double *, RESULT); // entkommentieren

  
#if defined  OBSOLETE_RFU
  //  #warning obsolete compilation
  #undef AttachMessageN
  typedef void (*setparameterfct) (int, int, SEXP, char[200], bool, int);
  typedef void (*getparameterfct) (SEXP, int, int);
  typedef void (*finalsetparameterfct) (int);
  typedef void (*deleteparameterfct) (int);

  #define MAXNLIST 7
  
  DECLARE1(void, getUtilsParam, utilsoption_type **, up)
  DECLARE1(void, getErrorString, errorstring_type, errorstring)
  DECLARE1(void, setErrorLoc, errorloc_type, errorloc)
  DECLARE1(void, relaxUnknownRFoption, bool, relax)
  DECLARE10(void, attachRFoptions, const char **, prefixlist, int, N, 
	   const char ***, all, int *, allN, setparameterfct, set, 
	   finalsetparameterfct, final, getparameterfct, get,
	    deleteparameterfct, del,
	   int, pl_offset,
	   bool, basicopt)
  DECLARE2(void, detachRFoptions, const char **, prefixlist, int, N)
  DECLARE3(int *, ToIntI, SEXP, X, bool *, create, bool, round)
  #define warn(X) {RFWARNING(X);}

#endif

  

 /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
  */
 
  /*

    See in R package RandomFields, /src/userinterfaces.cc 
          CALL#(...)
    at the beginning for how to make the functions available
    in a calling package

  */
#ifdef __cplusplus
}
#endif


/*

install.packages("RandomFieldsUtils_0.5.21.tar.gz", configure.args="CXX_FLAGS=-march=native", repo=NULL); library(RandomFieldsUtils, verbose=100); q()

*/

#define noMISS 0
#define noUSE 0
#define anyrelevantUSE 0
#define gpuUSE 1
#define avx2USE 2
#define avxUSE 3 
#define ssse3USE 4
#define sse2USE 5
#define avx512fUSE 6
#define USEnMISS 10
#define gpuMISS 11
#define avx2MISS 12
#define avxMISS 13 
#define ssse3MISS 14
#define sse2MISS 15
#define avx512fMISS 16
#define anyMISS (1 << gpuMISS) | (1 << avx2MISS) | (1 << avxMISS) | \
  (1 << ssse3MISS) | (1 << sse2MISS) | (1 << avx512fMISS)
#define AVX_INFO					\
  allmiss | alluse | (HAS_PARALLEL || alluse != 0) * (1 << anyrelevantUSE) | \
  ((HAS_PARALLEL || alluse != noUSE) && !(HAS_PARALLEL && allmiss==noMISS)) * \
  (1 << USEnMISS)



#define EAX 0
#define EBX 1
#define ECX 2
#define EDX 3
#define sse Available(1, EDX,25)
#define sse2 Available(1, EDX,26)
#define sse3 Available(1, ECX,0)
#define ssse3 Available(1, ECX,9)
#define sse41 Available(1, ECX,19)
#define avx Available(1, ECX,28)
#define avx2 Available(7, EBX,5)
#define avx512f Available(7, EBX,16)
#define avx512df Available(7, EBX,26)
#define avx512ar Available(7, EBX,27)
#define avx512cd Available(7, EBX,28)

#if defined WINCPUID
  #define AVAILABLE_SIMD static inline bool			\
       Available(unsigned Blatt, int Register, int Bit) {\
      uint32_t s[4];							\
      __cpuid((int *)s, (int) Blatt);					\
      return s[Register] & (1 << (Bit));				\
    }
  #if defined WIN32 || defined _WIN32 || defined __WIN32__
  #else
    #error Puzzled about the underlying system. Please contact maintainer.
  #endif
#elif defined LINUXCPUID
  #define AVAILABLE_SIMD static inline bool				\
    Available(unsigned Blatt, int Register, int Bit) {			\
      uint32_t s[4];							\
      asm volatile							\
      ("cpuid": "=a"(s[0]), "=b"(s[1]),"=c"(s[2]),			\
	"=d"(s[3]):"a"(Blatt),"c"(0));					\
      return s[Register] & (1 << (Bit));				\
      }
#else
  #define AVAILABLE_SIMD static inline bool				\
      Available(unsigned VARIABLE_IS_NOT_USED B, int VARIABLE_IS_NOT_USED R, \
		int VARIABLE_IS_NOT_USED Bit) {			\
         ERR("SIMD checks are not available on your system (on MS systems only under Visual Studio). Use 'SERVER' on Linux systems and alike.");  \
      return false;							\
      }						
#endif


#define WARN_PARALLELXX							\
  if (OPTIONS.basic.warn_parallel && mypid == parentpid) 	\
    PRINTF("Do not forget to run 'RFoptions(storing=FALSE)' after each call of a parallel command (e.g. from packages 'parallel') that calls a function in 'RandomFields'. (OMP within RandomFields is not affected.) This message can be suppressed by 'RFoptions(warn_parallel=FALSE)'.") /*// ok */ \
    
#define WARN_PARALLEL 

#define ASSERT_AVAILABILITY(V,W) if (!(V)) {char msg[300]; SPRINTF(msg, ASSERT_TEXT, #V, STRCMP(#V, #W) ? " && " : "", STRCMP(#V, #W) ? #W : ""); RFERROR(msg);}
#if defined SERVER_CAPACITY
//  #warning SERVER_CAPACITY 
  #if defined REQUIRED_SIMD
//    #warning SERVER = TRUE
   #define ASSERT_TEXT							\
"The program was compiled for '%.10s%.5s%.10s', but the CPU doesn't know about it. As 'SERVER=TRUE' has been chosen as compilation option, it was assumed that the programme was compiled on the most unskilled CPU." // ok
  #else
//  #warning SERVER = avx, sse, etc
    #define ASSERT_TEXT \
   "The program was compiled for '%.10s%.5s%.10s', but the CPU doesn't know about it. As 'SERVER' has been chosen as compilation option, it was assumed that each CPU has at least the SERVER skills."
   #endif
#elif defined REQUIRED_SIMD // ! SERVER_CAPACITY
//  #warning NO SERVER_CAPACITY given
  #if REQUIRED_SIMD == 0 // SERVER = nosimd without -mno-sse2
    #define ASSERT_TEXT \
    "The program was compiled for '%.10s%.5s%.10s' (what should mean 'without SIMD'), but the compiler includes SIMD. ('SERVER=nosimd' has been chosen.)"
  #elif REQUIRED_SIMD == 1 // SERVER = nosimd and -mno-sse2
    #define ASSERT_TEXT\
    "The program was compiled for '%.10s%.5s%.10s' (what should mean 'without SIMD', but the CPU requires SIMD at a higher level.  Please contact the maintainer."
  #elif REQUIRED_SIMD == 2 // SERVER=NA
    #define ASSERT_TEXT \
    "The program was compiled for '%.10s%.5s%.10s' (what should mean 'without SIMD'), but the compiler includes SIMD (at a higher level). ('SERVER=NA' had been chosen.)"
  #elif REQUIRED_SIMD == 3 // SERVER=noflags, FALSE
    #undef ASSERT_AVAILABILITY
    #define ASSERT_AVAILABILITY(V,W)
  #else
    #define ASSERT_TEXT\
    "The program was compiled for '%.10s%.5s%.10s' and leads to an unexpected situation! Please contact the maintainer."
  #endif
#else
  #define ASSERT_TEXT\
     "The program was compiled for '%.10s%.5s%.10s' and leads to an unexpected situation. Please contact the maintainer."
#endif



#if defined AVX512F
#define ASSERT_SIMD_CONSISTENCY(V) ASSERT_AVAILABILITY(avx512f, V)
#define ASSERT_SIMD_CONSISTENCY0 ASSERT_AVAILABILITY(avx512f, avx512f)
#elif defined AVX2
#define ASSERT_SIMD_CONSISTENCY(V) ASSERT_AVAILABILITY(avx2, V)
#define ASSERT_SIMD_CONSISTENCY0 ASSERT_AVAILABILITY(avx2, avx2)
#elif defined AVX
#define ASSERT_SIMD_CONSISTENCY(V) ASSERT_AVAILABILITY(avx, V)
#define ASSERT_SIMD_CONSISTENCY0 ASSERT_AVAILABILITY(avx, avx)
#elif defined SSSE3
#define ASSERT_SIMD_CONSISTENCY(V) ASSERT_AVAILABILITY(ssse3, V)
#define ASSERT_SIMD_CONSISTENCY0 ASSERT_AVAILABILITY(ssse3, ssse3)
#elif defined SSE3
#define ASSERT_SIMD_CONSISTENCY(V) ASSERT_AVAILABILITY(sse3, V)
#define ASSERT_SIMD_CONSISTENCY0 ASSERT_AVAILABILITY(sse3, sse3)
#elif defined SSE2
#define ASSERT_SIMD_CONSISTENCY(V) ASSERT_AVAILABILITY(sse2, V)
#define ASSERT_SIMD_CONSISTENCY0 ASSERT_AVAILABILITY(sse2, sse2)
#else
#define ASSERT_SIMD_CONSISTENCY(V) ASSERT_AVAILABILITY(V, V)
#define ASSERT_SIMD_CONSISTENCY0 
#endif
  

#define C_ASSERT_ANYSIMD(FILE)				\
    Uint alluse = noUSE,				\
      miss = noMISS,			\
      allmiss = noMISS;			\
    check_simd_##FILE()

#define C_ASSERT_SIMD(FILE, WHAT)		\
    allmiss |= (miss = check_simd_##FILE());\
    if (!miss) alluse |= 1 << WHAT##USE;
    
    
#define D_ASSERT_ANYSIMD(FILE) extern Uint check_simd_##FILE() // ok dass nicht in extern.h:
//#define D_ASSERT_SIMD(FILE, WHAT) D_ASSERT_ANYSIMD(FILE)
#define D_ASSERT_SIMD(FILE) D_ASSERT_ANYSIMD(FILE)

#define ASSERT_ANYSIMD(FILE) Uint check_simd_##FILE() {ASSERT_SIMD_CONSISTENCY0; return noMISS;}
#define ASSERT_SIMD(FILE, WHAT) Uint check_simd_##FILE() { ASSERT_SIMD_CONSISTENCY(WHAT); return noMISS;}

#define ASSERT_SIMD_MISS(FILE, WHAT) Uint check_simd_##FILE() { return 1<<WHAT##MISS; }

#endif
