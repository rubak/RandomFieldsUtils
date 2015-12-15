


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 Martin Schlather

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



#ifndef rfutils_solve_H
#define rfutils_solve_H 1


typedef enum InversionMethod { Cholesky, SVD, Sparse, QR, LU, 
			       NoFurtherInversionMethod, NoInversionMethod,
			       Diagonal // always last one!
} 
InversionMethod;
extern const char * InversionNames[(int) Diagonal + 1];

#define PIVOT_NONE 0
#define PIVOT_MMD 1
#define PIVOT_RCM 2

#define SOLVE_SVD_TOL 3
#define SOLVE_METHODS 3
#define solveN 11
typedef struct solve_param {
  double sparse, spam_tol, spam_min_p, svd_tol;
  InversionMethod Methods[SOLVE_METHODS];
  int spam_min_n, spam_sample_n, spam_factor,
    pivot, max_chol, max_svd;
} solve_param;

#define SOLVE_NAMES \
  { "use_spam", "spam_tol", "spam_min_p", "svdtol",			\
    "matrix_methods", "spam_min_n", "spam_sample_n", "spam_factor",	\
      "spam_pivot", "max_chol", "max_svd"				\
  }


#define solve_param_default				\
  { RF_NA, DBL_EPSILON,	0.8, 1e-8,			\
      {NoInversionMethod, NoInversionMethod},		\
      400, 500, 4294967, PIVOT_MMD, 1000000000, 1000000000}


#define solve_param_default_RF				\
  { RF_NA, DBL_EPSILON,	0.8, 1e-8,			\
      {NoInversionMethod, NoInversionMethod},		\
      400, 500, 4294967, PIVOT_MMD, 8192, 6555}

extern solve_param SolveParam;

typedef struct solve_storage {
  int SICH_n, MM_n, workspaceD_n, workspaceU_n, VT_n, U_n, D_n, 
    iwork_n, work_n, w2_n, ipiv_n, workLU_n, pivot_n,
    xlnz_n, snode_n, xsuper_n, xlindx_n, invp_n, 
    cols_n, rows_n, DD_n, lindx_n, xja_n,
    lnz_n, RHS_n, w3_n;
  //t_cols_n, t_rows_n, t_DD_n;
  InversionMethod method, newMethods[SOLVE_METHODS];
  int nsuper, nnzlindx, size,
    *iwork, *ipiv,
    *pivot, *xlnz, *snode, *xsuper, *xlindx, 
    *invp, *cols, *rows, *lindx, *xja; //*t_cols, *t_rows;
  double *SICH, *MM, *workspaceD, *workspaceU,
    *VT, *work, *w2, *U, *D, *workLU, 
    *lnz, *DD, *w3, *RHS; //, *t_DD;
} solve_storage;



#define CMALLOC(WHICH, N, TYPE)	{				 \
    int _N_ = N;						 \
    if (pt->WHICH##_n < _N_) {					 \
      if (pt->WHICH##_n < 0) BUG;				 \
      FREE(pt->WHICH);						 \
      pt->WHICH##_n = _N_;						\
	if ((pt->WHICH = (TYPE *) CALLOC(_N_, sizeof(TYPE))) == NULL)	\
	  return ERRORMEMORYALLOCATION;					\
    } else {								\
      assert( (_N_ > 0 && pt->WHICH != NULL) || _N_ == 0);		\
      for (int iii=0; iii<_N_; pt->WHICH[iii++] = 0);			\
    }									\
  }									\
    TYPE VARIABLE_IS_NOT_USED *WHICH = pt->WHICH			



#endif
