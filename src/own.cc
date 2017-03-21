/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- Martin Schlather, Reinhard Furrer

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

//#include "Basic_utils.h"  // must be before anything else
#include "RandomFieldsUtils.h"  // must be before anything else
#ifdef DO_PARALLEL
#include <omp.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include "General_utils.h"
#include "own.h"
#include "init_RandomFieldsUtils.h"



// local
char ERRMSG[LENERRMSG], MSG[LENERRMSG], BUG_MSG[250], MSG2[LENERRMSG];

// globally needed
errorloc_type ERROR_LOC="";
errorstring_type ERRORSTRING;


SEXP attachRFoptionsUtils() {
  //  NList = 0;
  
  // printf("UTx %ld\n", (long) getUtilsParam);

  attachRFoptions(ownprefixlist, ownprefixN, ownall, ownallN,
  		  setparameterUtils, NULL, getparameterUtils);

#ifdef DO_PARALLEL
  basic_param *gp = &(GLOBAL.basic);
  omp_set_num_threads(gp->cores);
#endif
  
  return R_NilValue;
}

SEXP detachRFoptionsUtils(){
#ifdef DO_PARALLEL      
  omp_set_num_threads(1);
#endif
  detachRFoptions(ownprefixlist, ownprefixN);
  return R_NilValue;
}
