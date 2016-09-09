
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 Martin Schlather

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

#include "Basic_utils.h" // must be before anything else
#ifdef DO_PARALLEL
#include <omp.h>
#endif
#include "General_utils.h"
#include "kleinkram.h"
#include "init_RandomFieldsUtils.h"

// IMPORTANT: all names of general must be at least 3 letters long !!!
const char *basic[basicN] =
  { "printlevel", "skipchecks", "cPrintlevel", "seed",  "asList", "cores"};

const char * solve[solveN] = 
  { "use_spam", "spam_tol", "spam_min_p", "svdtol",			
    "solve_method", "spam_min_n", "spam_sample_n", "spam_factor",	
    "spam_pivot", "max_chol", "max_svd", "eigen2zero"
    //, "tmp_delete"
  };


#define ownprefixN 2
const char * ownprefixlist[ownprefixN] = {"basic", "solve"};
const char **ownall[ownprefixN] = {basic, solve};
int ownallN[ownprefixN] = {basicN, solveN};


int PL=C_PRINTLEVEL;

utilsparam GLOBAL = {
  basic_START,
  solve_START
};







void setparameterUtils(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		       bool VARIABLE_IS_NOT_USED isList) {
  switch(i) {
  case 0: {// general
    basic_param *gp;
    gp = &(GLOBAL.basic);
    switch(j) {
    case 0: {
      int threshold = 1000; //PL_ERRORS;
      gp->Rprintlevel = INT;
      PL = gp->Cprintlevel = 
	gp->Rprintlevel <= threshold ? gp->Rprintlevel : threshold;
    }
      break;
    case 1: gp->skipchecks = LOG;    break;
    case 2: PL = gp->Cprintlevel = INT;        break;
    case 3: gp->seed = Integer(el, name, 0, true); break;
    case 4: gp->asList = LOG; break;
    case 5: gp->cores = POSINT; 
#ifdef DO_PARALLEL      
      omp_set_num_threads(gp->cores);
#else
      if (gp->cores != 1) 
	ERR("The system does not allow for OpenMP: the value 1 is kept for'cores'.");
#endif
      break;
    default: BUG;
    }}
    break;
 
 case 1: {
    solve_param *so = &(GLOBAL.solve);
    switch(j) {
    case 0: {
      double sparse = NUM;
      so->sparse = !R_finite(sparse) ? Nan : sparse==0.0 ? False : True ; 
      break;
    }
    case 1: so->spam_tol = POS0NUM; break;      
    case 2: so->spam_min_p = POS0NUM; break;      
    case SOLVE_SVD_TOL: so->svd_tol = POS0NUM; break;        
    case 4: GetName(el, name, InversionNames, nr_user_InversionMethods,
		    (int) NoInversionMethod, (int) NoFurtherInversionMethod, 
		    (int *)so->Methods, SOLVE_METHODS);
      break;
    case 5: so->spam_min_n = POSINT; break;      
    case 6: so->spam_sample_n = POSINT; break;      
    case 7: so->spam_factor = POSINT; break;      
    case 8: so->pivot = POSINT; 
      if (so->pivot > 2) so->pivot = PIVOT_NONE;
      break;    
    case 9: so->max_chol = POSINT; break;      
    case 10: so->max_svd = POSINT; break;    
      //    case 11: so->tmp_delete = LOG; break;    
    case 11: so->eigen2zero = POS0NUM; break;        
    default: BUG;
    }}
    break;
    
  default: BUG;
  }

}


void getparameterUtils(SEXP *sublist) {
  int i, k;
    //#define ADD(ELT) {printf(#ELT"\n");SET_VECTOR_ELT(sublist[i], k++, ELT);}
    i = 0; { 
    // printf("OK %d\n", i);
    k = 0;
    basic_param *p = &(GLOBAL.basic);
    ADD(ScalarInteger(p->Rprintlevel));    
    ADD(ScalarLogical(p->skipchecks));
    ADD(ScalarInteger(p->Cprintlevel));
    ADD(ScalarInteger(p->seed));    
    ADD(ScalarLogical(p->asList));   
    ADD(ScalarInteger(p->cores));    
   }
  
 i++; {
    k = 0;
    solve_param *p = &(GLOBAL.solve);
    // printf("sparse user interface %d %d; %d %d\n", p->sparse, ExtendedBoolean(p->sparse), NA_LOGICAL, NA_INTEGER);
    ADD(ExtendedBooleanUsr(p->sparse));
    ADD(ScalarReal(p->spam_tol));    
    ADD(ScalarReal(p->spam_min_p));    
    ADD(ScalarReal(p->svd_tol)); 
    SET_VECTOR_ELT(sublist[i], k++,
		   String((int*) p->Methods, InversionNames, SOLVE_METHODS,
			  (int) NoFurtherInversionMethod));	
    ADD(ScalarInteger(p->spam_min_n));    
    ADD(ScalarInteger(p->spam_sample_n));    
    ADD(ScalarInteger(p->spam_factor));    
    ADD(ScalarInteger(p->pivot));    
    ADD(ScalarInteger(p->max_chol));    
    ADD(ScalarInteger(p->max_svd)); 
    ADD(ScalarReal(p->eigen2zero)); 
    //    ADD(ScalarLogical(p->tmp_delete));
  }
 
 
 assert (i == ownprefixN - 1); 
}


