/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 -- 2021 Martin Schlather

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

#include <unistd.h>
#include "parallel_simd.h"
#ifdef TIME_AVAILABLE
#include <time.h>
#endif

#include "Basic_utils_local.h" // must be before anything else
#include "kleinkram.h"
#include "zzz_RandomFieldsUtils.h"
#include "xport_import.h"
#include "options.h"
#include "extern.h"


#define PLverbose 2

// IMPORTANT: all names of general must have at least 3 letters !!!
const char *basic[basicN] =
  { "printlevel","cPrintlevel", "seed", "cores",
    "skipchecks", "asList", "verbose", "helpinfo", "efficient"
  };

const char *installNrun[installNrunN] =
  { "kahanCorrection", "warn_unknown_option", "la_mode", "warn_parallel",
    "install","installPackages", "determineLAmode", "mem_is_aligned",
    "gpuDevices", "maxStreams"
  };

const char * solve[solveN] = 
  { "use_spam", "spam_tol", "spam_min_p", "svdtol", "eigen2zero",
    "solve_method", "spam_min_n", "spam_sample_n", "spam_factor", "spam_pivot",
    "max_chol", "max_svd", "pivot",
    "pivot_idx", // dynamic parameter
    "pivot_relerror", "pivot_max_deviation", "pivot_max_reldeviation",
    "det_as_log", "pivot_actual_size", "pivot_check", "pseudoinverse"
    //, "tmp_delete"
  };


const char * prefixlist[prefixN] = {"basic", "installNrun", "solve"};
const char **allOptions[prefixN] = {basic, installNrun, solve};
int allOptionsN[prefixN] = {basicN, installNrunN, solveN};


utilsoption_type OPTIONS = { // OK
  basic_START,
  installNrun_START,
  { solve_START }
};


//#if defined(unix) || defined(__unix__) || defined(__unix)
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
#else
int numCPU = MAXINT;
#endif



int doPosDefIntern(double *M0, int size, bool posdef,
		   double *rhs0, Long rhs_cols, double *result,
		   double *logdet, int calculate, solve_storage *Pt,
		   solve_options *sp, int VARIABLE_IS_NOT_USED cores);



void setoptionsRFU(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		   bool isList, utilsoption_type *options) {
  switch(i) {
  case 0: {// general
    basic_options *gp = &(options->basic);
    switch(j) {
    case 0: {  // general options
      int threshold = 1000; //PL_ERRORS;
      gp->Rprintlevel = INT;
      if (gp->Rprintlevel > threshold) gp->Rprintlevel = threshold;
      PL = gp->Cprintlevel = gp->Rprintlevel + PLoffset;
    }
      break;
    case 1: PL = gp->Cprintlevel = INT + PLoffset ;
      break;
    case 2: gp->seed = Integer(el, name, 0, true); break;
    case 3: gp->cores = POSINT;
      if (gp->cores > numCPU) {
	WARN1("number of 'cores' is set to %d", numCPU);
	gp->cores = numCPU;
      }
#ifdef DO_PARALLEL
#else
      if (gp->cores != 1) {
	gp->cores = 1;
	PRINTF("The system does not allow for OpenMP: the value 1 is kept for 'cores'.");
      }    
#endif
      if (options->installNrun.determineLAmode) SetLaMode();
      CORES = gp->cores;
      break;
    case 4: gp->skipchecks = LOGI;    break;
    case 5: gp->asList = LOGI; break;
    case 6 : if (!isList) {
	PL = gp->Cprintlevel = gp->Rprintlevel = 1 + (LOGI * (PLverbose - 1));
      }
      break;
    case 7: gp->helpinfo = LOGI; break;
    case 8 : gp->efficient = LOGI; break;
    default: BUG;
    }}
    break;
    
  case 1: {
    installNrun_options *gp = &(options->installNrun);
    switch(j) {
    case 0: gp->kahanCorrection = LOGI; break;
    case INSTALL_RUN_WARN_OPTION: gp->warn_unknown_option = INT;
      break;
    case 2: {
      int neu;
      if (TYPEOF(el) == STRSXP) {
	neu = GetName(el, name, LA_NAMES, LA_LAST + 1, gp->la_usr);
	//if (neu == LA_QUERY) {
	// if (!local) // OK
	//   PRINTF("internal mode = '%.10s'\nmax size for internal Cholesky = %d\nmax size for tiny size algorithm = %d\n", LA_NAMES[gp->la_mode], gp->LaMaxTakeIntern,  options->solve.tinysize);
	//}
      } else {
	neu = POS0INT;
	if (neu > LA_LAST) ERR0("wrong value for 'la_mode'");
      }				  
      if (neu != LA_QUERY) {
#ifdef USEGPU
#else
	if (neu > LA_GPU)
	  ERR1("In order to use '%.20s', install the package with apropriate compiler options.",
	       LA_NAMES[LA_GPU]);
#endif
	SetLaMode((la_modes) neu,
		  options->basic.cores);  // koennen noch fehler auftreten
	//      printf("simu ende\n");
	gp->la_usr = (la_modes) neu;
      }
    }
      break;     
    case 3 : gp->warn_parallel = LOGI;  break;
    case 4 : {
      install_modes old =  gp->install;
      gp->install = (install_modes )
	GetName(el, name, INSTALL_NAMES, INSTALL_LAST + 1, Iask);
      if (gp->install == Inone) gp->installPackages = false;
      else if (old != gp->install) {
	gp->installPackages = true;
	resetInstalled();
      }
    }
      break;
    case 5 :
      // gp->installPackages != LOGI;  
      break;
    case 6 : gp->determineLAmode = LOGI;  break;
    case 7 :
      // gp->mem_is_aligned = LOGI;
      break;
    case 8 :  Integer(el, name, gp->gpu_devices, MAX_GPU_DEVICES) ;
      gp-> Ngpu_devices = MIN(length(el), MAX_GPU_DEVICES);
      break;
     case 9 :
       gp->maxStreams = POS0INT;
      break;
   default: BUG;
    }}
    break;
    
  case 2: {
    //   printf("name = %.50s %d\n", name, j);
    
    solve_options *so = &(options->solve);
    switch(j) {
    case 0: so->sparse = USRLOG;
      if (so->sparse != False) {
	so->sparse = False;
	ERR0("'spam' is currently disabled.")
	  }      
      break; // USRLOGRELAXED??
    case 1: so->spam_tol = POS0NUM; break;
    case 2: Real(el, name, so->spam_min_p, 2);
      for (int u=0; u<=1; u++)
	so->spam_min_p[u] = so->spam_min_p[u] < 0 ? 0
  	: so->spam_min_p[u] > 1.0 ? 1.0 : so->spam_min_p[u];
      break;      
      case SOLVE_SVD_TOL: so->svd_tol = POS0NUM; break;        
    case 4: so->eigen2zero = POS0NUM; break;        
    case 5:
      GetName(el, name, InversionNames, nr_user_InversionMethods,
		    (int) NoInversionMethod, (int) NoFurtherInversionMethod, 
		    (int *)so->Methods, SOLVE_METHODS);
      break;
    case 6: Integer(el, name, so->spam_min_n, 2);      
    break; 
    case 7: so->spam_sample_n = POSINT; break;      
    case 8: so->spam_factor = POSINT; break;      
    case 9: so->pivotsparse = POSINT; 
      if (so->pivotsparse > 2) so->pivotsparse = PIVOT_NONE;
      break;    
    case 10:
      //      printf("max chol = %d\n",  so->max_chol);
      so->max_chol = POSINT;
      //      printf("X max chol = %d\n",  so->max_chol);
       break;      
    case 11: so->max_svd = POS0INT; break;    
      //    case 11: so->tmp_delete = LOGI; break;    
    case 12: 
      so->pivot_mode = (pivot_modes) GetName(el, name, PIVOT_NAMES,
					     PIVOT_LAST + 1, so->pivot_mode);
      break;    
    case 13: if (!isList) {
      int len = length(el);
      if (len == 0) {
	if (so->pivot_idx_n > 0) { FREE(so->pivot_idx); }
      } else {
	if (so->pivot_idx_n != len) {
	  FREE(so->pivot_idx);
	  so->pivot_idx = (int*) MALLOC(len * sizeof(int));
	}
	for (int L=0; L<len; L++) so->pivot_idx[L] = Integer(el, name, L);
      }
      so->pivot_idx_n = len;     
    }
      break; 
    case 14: so->pivot_relerror = POS0NUM; break;    
    case 15: so->max_deviation = POSNUM; break;    
    case 16: so->max_reldeviation = POS0NUM; break;    
    case 17: so->det_as_log = LOGI; break;    
    case 18: so->actual_size = POS0NUM; break;    
    case 19: so->pivot_check = USRLOG; break;    
    case 20: so->pseudoinverse = LOGI; break;
    default: BUG;
    }}
    break;
    
    default: BUG;
  }
}

void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		       bool isList, bool local) {
  if (!local && parallel()) 
    ERR1("Option '%.25s' can be set only through 'RFoptions' at global level",
	 allOptions[i][j]);
  setoptionsRFU(i, j, el, name, isList,  WhichOptionList(local));
}



void getoptionsRFU(SEXP sublist, int i, utilsoption_type *options) {  
   int  k = 0;
 //printf("OK %d\n", i);    
 switch(i) {
  case 0 : {
    //    printf("OK %d\n", i);
    basic_options *p = &(options->basic);
    ADD(ScalarInteger(p->Rprintlevel));    
    ADD(ScalarInteger(p->Cprintlevel - PLoffset));
    ADD(ScalarInteger(p->seed));    
    ADD(ScalarInteger(p->cores));    
    ADD(ScalarLogical(p->skipchecks));
    ADD(ScalarLogical(p->asList));   
    ADD(ScalarLogical(p->Rprintlevel >= PLverbose));
    ADD(ScalarLogical(p->helpinfo));    
    ADD(ScalarLogical(p->efficient));
   }
    break;
 
  case 1 : {
    installNrun_options *p = &(options->installNrun);
    ADD(ScalarLogical(p->kahanCorrection));   
    ADD(ScalarInteger(p->warn_unknown_option));    
    ADD(ScalarString(mkChar(LA_NAMES[p->la_usr])));
    ADD(ScalarLogical(p->warn_parallel));
    ADD(ScalarString(mkChar(INSTALL_NAMES[p->install])));
    ADD(ScalarLogical(p->installPackages));
    ADD(ScalarLogical(p->determineLAmode));
    ADD(ScalarLogical(p->mem_is_aligned));
    SET_VECTOR_ELT(sublist, k++, Int(p->gpu_devices, p->Ngpu_devices)); 
    ADD(ScalarInteger(p->maxStreams));    
  }
    break;
 
  case 2 : {
    solve_options *p = &(options->solve);
    //    printf("sparse user interface %d %d %d\n", p->sparse, NA_LOGICAL, NA_INTEGER);
    ADD(ExtendedBooleanUsr(p->sparse));
    //
      ADD(ScalarReal(p->spam_tol));    
    SET_VECTOR_ELT(sublist, k++, Num(p->spam_min_p, 2));
   ADD(ScalarReal(p->svd_tol)); 
    ADD(ScalarReal(p->eigen2zero));
      
   SET_VECTOR_ELT(sublist, k++,
		   String((int*) p->Methods, InversionNames, SOLVE_METHODS,
			  (int) NoFurtherInversionMethod));	
   //   printf("A\n");
   SET_VECTOR_ELT(sublist, k++, Int(p->spam_min_n, 2));
   ADD(ScalarInteger(p->spam_sample_n));    
    ADD(ScalarInteger(p->spam_factor));    
    ADD(ScalarInteger(p->pivotsparse));    
    ADD(ScalarInteger(p->max_chol));    
    ADD(ScalarInteger(p->max_svd)); 
    ADD(ScalarString(mkChar(PIVOT_NAMES[p->pivot_mode])));
    //if (true)
    SET_VECTOR_ELT(sublist, k++, Int(p->pivot_idx, p->pivot_idx_n));
    //  else ADD(ScalarInteger(NA_INTEGER));
    //    ADD(ScalarLogical(p->tmp_delete));
     ADD(ScalarReal(p->pivot_relerror));    
    ADD(ScalarReal(p->max_deviation));    
    ADD(ScalarReal(p->max_reldeviation));    
    ADD(ScalarLogical(p->det_as_log));
    ADD(ScalarInteger(p->actual_size));
    ADD(ExtendedBooleanUsr(p->pivot_check));
    ADD(ScalarLogical(p->pseudoinverse));
  }
    break;
  default : BUG;
  }
  //  printf("EE A\n");
}

void getoptions(SEXP sublist, int i, bool local) {  
  getoptionsRFU(sublist, i, WhichOptionList(local));
}

void params_utilsoption(int local, int *params) {
  utilsoption_type *from = &OPTIONS;
  if (local) {
    KEY_type *KT = KEYT();
    from = &(KT->global_utils);
  } 
  params[PIVOT_IDX_N] = from->solve.pivot_idx_n;
}

void get_utilsoption(utilsoption_type *S, int local) {
  assert(solveN == 21 && basicN == 9 && installNrunN == 10 && prefixN==3);
  utilsoption_type *from = &OPTIONS;
  if (local) {
    KEY_type *KT = KEYT();
    from = &(KT->global_utils);
  }
  assert(from->solve.pivot_idx_n!=0 xor from->solve.pivot_idx == NULL);
  assert(S->solve.pivot_idx_n!=0 xor S->solve.pivot_idx == NULL);
  if (S->solve.pivot_idx_n != from->solve.pivot_idx_n) BUG;
  int *save_idx = S->solve.pivot_idx;
  MEMCOPY(S, from, sizeof(utilsoption_type)); // OK
  S->solve.pivot_idx = save_idx;
  if (S->solve.pivot_idx_n > 0) {
    MEMCOPY(S->solve.pivot_idx, from->solve.pivot_idx, 
	    sizeof(int) * S->solve.pivot_idx_n);
  }
}

void push_utilsoption(utilsoption_type *S, int local) {
  utilsoption_type *to = &OPTIONS;
  if (local) {
    KEY_type *KT = KEYT();
    to = &(KT->global_utils);
  } 
  assert(to->solve.pivot_idx_n!=0 xor to->solve.pivot_idx == NULL);
  assert(S->solve.pivot_idx_n!=0 xor S->solve.pivot_idx == NULL);
  int *save_idx = to->solve.pivot_idx;
  if (to->solve.pivot_idx_n != S->solve.pivot_idx_n) {
    FREE(to->solve.pivot_idx);
    to->solve.pivot_idx = (int*) MALLOC(S->solve.pivot_idx_n * sizeof(int));
    save_idx = to->solve.pivot_idx;
  }
  MEMCOPY(to, S, sizeof(utilsoption_type)); // OK
  to->solve.pivot_idx = save_idx;
  if (S->solve.pivot_idx_n > 0) {
    MEMCOPY(to->solve.pivot_idx, S->solve.pivot_idx,
	    sizeof(int) * S->solve.pivot_idx_n);
  }
}

void del_utilsoption(utilsoption_type *S) {
  FREE(S->solve.pivot_idx);
  S->solve.pivot_idx_n = 0;
}



  extern bool obsolete_package_in_use;

#define FASTER 1.3 // 1.3 as for typical application of likelihood,
// the determinant calculation in RandomFieldsUtils is for free. Somehow a
// balance
int own_chol_up_to(int size, int maxtime, int VARIABLE_IS_NOT_USED cores) {
#ifdef TIME_AVAILABLE
  if (size <= 0) return true;
  Long delta[2];
  solve_options sp;
  solve_storage pt;
  solve_NULL(&pt);
  MEMCOPY(&sp, &(OPTIONS.solve), sizeof(solve_options)); 
  sp.Methods[0] = Cholesky;					     
  sp.Methods[1] = NoFurtherInversionMethod;
  sp.pivot_mode = PIVOT_NONE;	  
  sp.sparse = False;	  
  double old_quotient = RF_NAN;
  // basic assumption is that R implementation's getting faster
  // for larger matrices
  //  printf("**** start\n");
  while (true) {
    //    printf("x \n");
    int sizeP1 = size + 1,
      sizesq = size * size,
      loops = size > 64 ? 1 : 16384 / ((size + 8) * (size+8)) / 4;
    double *M = (double*) MALLOC(sizesq * sizeof(double));
    for (int j=0; j<=1; j++) {
      //      printf("%d,%d\n", j, loops);
      SetLaMode(j == 0
		|| obsolete_package_in_use
		? LA_INTERN : LA_R, cores);
      clock_t start = clock();
      for (int k=0; k<loops; k++) {      
	MEMSET(M, 0, sizesq * sizeof(double));
	for (int i=0; i<sizesq; i+=sizeP1) M[i] = 1.0;
	if (size > 1) M[1] = M[size] = 1e-5;
	//printf("size=%d\n", size);
	doPosDefIntern(M, size, true, NULL, 0, NULL, NULL, MATRIXSQRT, &pt, &sp,
		       cores);
	//printf("doen\n");
      }
      delta[j] = (Long) clock() - start; 
      if (delta[j] < 0) delta[j] += MAXINT; // manual: 32bit repres. of clock
    }
    FREE(M);
    if (PL > 2)
      PRINTF("Cholesky decomposition for a %d x %d matrix needs %ld and %ld [mu s] on R and facile code on %d cores (#%d), respectively.\n", size, size, delta[1], delta[0], CORES, loops);

    //  printf("delta %d %d %d\n", delta[0], delta[1], maxtime);
    
    if (delta[0] > maxtime || delta[1] > maxtime ||
	delta[0] >= FASTER * delta[1]){
      solve_DELETE0(&pt);
      if (maxtime > 0 &&
	  (delta[0] > 10 * maxtime || delta[1] > 10 * maxtime) ||
	  delta[0] > 2 * delta[1] || delta[1] > 2 * delta[0]) {
	// seems to be time consuming. So stop.
	return (double) delta[0] < FASTER * (double) delta[1]
	  ? MAXINT : (size <= 0 ? 0 : size / 2);
      }
      break;
    }
    old_quotient = (double) delta[0] / delta[1];
    size *= 2;
  }
  double new_quotient = (double) delta[0] / delta[1];
  if (new_quotient < FASTER) return MAXINT;
  if (size <= 0) return(0);
  if (!R_FINITE(old_quotient)) {
    // printf("halfsize\n");
    int compare = own_chol_up_to(size / 2, 0, cores);
   return compare == MAXINT ? size : compare;
  }
  double x0 = 0.5 * size * (1.0 + (FASTER - old_quotient) /
			    (new_quotient - old_quotient)); //lin interpol
  assert(x0 >= 0.5 * size && x0 <= size);
  int compare = own_chol_up_to((int) x0, 0, cores);
  //  printf("%f %f %f %f %d %d\n", x0,FASTER,  old_quotient, new_quotient, size, compare);
  return (int) (compare == MAXINT ?  x0 : 0.5 * size);
#else 
  ERR0("option 'LA_AUTO' is available only on linux systems");
  return 0;
#endif
}
int own_chol_up_to(int VARIABLE_IS_NOT_USED cores) {
  own_chol_up_to(256, 0, cores); //warm up for some BLAS implementatioan
  //  CORES = GL OBAL.basic.cores = 4;
  //  own_chol_up_to(256, 50000);
  //  own_chol_up_to(8, 50000);
  return own_chol_up_to(256, 50000, cores);
}


  
void SetLaMode(la_modes usr_mode, int VARIABLE_IS_NOT_USED cores) {
  utilsoption_type *utils = &OPTIONS;
  la_modes la_mode = usr_mode;
  utils->solve.tinysize =
    utils->installNrun.LaMaxTakeIntern = -1;
#define TINY_SIZE_MAX 3  
  if (la_mode == LA_INTERN) {
    utils->solve.tinysize = TINY_SIZE_MAX;
    utils->installNrun.LaMaxTakeIntern = MAXINT;
  } else if (la_mode == LA_AUTO) {  
    la_mode = HAS_GPU ? LA_GPU : LA_R ;
#if defined TIME_AVAILABLE
#  ifdef SCHLATHERS_MACHINE
#else    
    int PLalt = PL;
    PL = 0;
#  endif  
    utils->installNrun.LaMaxTakeIntern = own_chol_up_to(cores);
    utils->solve.tinysize = MIN(TINY_SIZE_MAX, utils->installNrun.LaMaxTakeIntern);
    if (PL > 0)
      PRINTF("Limit size for facile Cholesky algorithm  = %d\n",
	     utils->installNrun.LaMaxTakeIntern);
#  ifdef SCHLATHERS_MACHINE
#else   
    PL = PLalt;
#  endif
#endif
  }
  
  if ((la_mode == LA_GPU || la_mode == LA_R) &&
      utils->solve.pivot_mode > PIVOT_AUTO)
    ERR0("Pivotized Cholesky decomposition has not been implemented yet for GPU and the LAPACK library");
  
  utils->installNrun.la_mode = la_mode;
}
void SetLaMode() {
  utilsoption_type *utils = &OPTIONS;
  SetLaMode(utils->installNrun.la_usr,
	    utils->basic.cores);
}


