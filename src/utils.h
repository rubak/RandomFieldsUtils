#ifndef rfutils_H
#define rfutils_H 1

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include "Basic_utils.h"
#include "Error_utils.h"
#include "Matrix_mult.h"
#include "Solve.h"
#include "Struve.h"


#define DOPRINTF if (DOPRINT) Rprintf
#define PRINTF Rprintf
#define print PRINTF /* // */


#ifndef SCHLATHERS_MACHINE
#define INTERNAL SERR("Sorry. This functionality does not exist currently. There is work in progress at the moment by the maintainer.")
#define assert(X) {}
#define VARIABLE_IS_NOT_USED
#define BUG {								\
    sprintf(BUG_MSG, "Severe error occured in function '%s' (file '%s', line %d). Please contact maintainer martin.schlather@math.uni-mannheim.de .", \
	    __FUNCTION__, __FILE__, __LINE__);				\
    error(BUG_MSG);							\
  }									
#define DO_TESTS false
#define ERRLINE 
#define MEMCOPY(A,B,C) memcpy(A,B,C)
#define MALLOC malloc
#define CALLOC calloc
#define FREE(X) if ((X) != NULL) {free(X); (X)=NULL;}
#define UNCONDFREE(X) {free(X); (X)=NULL;}


#else // SCHLATHERS_MACHINE

// __extension__ unterdrueckt Fehlermeldung wegen geklammerter Argumente
#define INTERNAL  \
  sprintf(BUG_MSG, \
	  "made to be an internal function '%s' ('%s', line %d).", /* // */ \
	  __FUNCTION__, __FILE__, __LINE__);				\
  warning(BUG_MSG)
 
#define assert(X) if (!__extension__ (X )) {				\
    sprintf(BUG_MSG,							\
	    "'assert(%s)' failed in function '%s' (file '%s', line %d).", \
	    #X,__FUNCTION__, __FILE__, __LINE__);			\
    error(BUG_MSG);							\
  }
#define VARIABLE_IS_NOT_USED __attribute__ ((unused))
#define SHOW_ADDRESSES 1
#define BUG {								\
    sprintf(BUG_MSG, "BUG in '%s' ('%s', line %d).", \
	    __FUNCTION__, __FILE__, __LINE__);				\
    error(BUG_MSG);							\
  }									
#define DO_TESTS true

#define MEMCOPY(A,B,C) ({ assert((A)!=NULL && (B)!=NULL); memcpy(A,B,C); })
//#define MEMCOPY(A,B,C) memory_copy(A, B, C)
#define MALLOC(X) ({assert(X>0 && X<=668467200);malloc(X);})
#define CALLOC(X, Y) ({assert((X)>0 && X<1e8 && (Y)>0 && (Y)<=64); calloc(X,Y);})
#define FREE(X) { if (showfree) DOPRINTF("(free in %s, line %d)\n", __FILE__, __LINE__); if ((X) != NULL) free(X); (X)=NULL;}
#define UNCONDFREE(X) { if (showfree) DOPRINTF("(free in %s, line %d)\n", __FILE__, __LINE__); free(X); (X)=NULL;}

#endif // SCHLATHERS_MACHINE



#ifdef RANDOMFIELDS_DEBUGGING

#undef MALLOC
#define MALLOC(X) ({DOPRINTF("(MALLOC %s, line %d)\n", __FILE__, __LINE__);assert(X>0 && X<=3e9);malloc(X);})
//
#undef CALLOC
#define CALLOC(X, Y) ({DOPRINTF("(CALLOC %s, line %d)\n",__FILE__, __LINE__);assert((X)>0 && X<1e8 && (Y)>0 && (Y)<=64); calloc(X,Y);})
//#define MALLOC malloc
//#define CALLOC calloc

#define DEBUGINFOERR {							\
    char dummy[MAXERRORSTRING]; strcpy(dummy, ERRORSTRING);		\
    sprintf(ERRORSTRING, "%s (%s, line %d)\n", dummy, __FILE__, __LINE__); \
  }
#define DEBUGINFO DOPRINTF("(currently at  %s, line %d)\n", __FILE__, __LINE__)

#else
#define DEBUGINFO 
#define DEBUGINFOERR
#endif


#define RF_NA NA_REAL 
#define RF_NAN R_NaN
#define RF_NEGINF R_NegInf
#define RF_INF R_PosInf
#define T_PI M_2_PI



#endif


