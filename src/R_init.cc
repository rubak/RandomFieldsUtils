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

#include <R.h>
#include <Rinternals.h>
#include "General_utils.h"
#include "R_init.h"
#include <R_ext/Rdynload.h>
#include "RandomFieldsUtils.h"


// local
char ERRMSG[LENERRMSG], MSG[LENERRMSG], BUG_MSG[250], MSG2[LENERRMSG];

// globally needed
char 
  ERRORSTRING[MAXERRORSTRING], 
  ERROR_LOC[nErrorLoc]="";


void getErrorString(char errorstring[MAXERRORSTRING]){
  strcopyN(errorstring, ERRORSTRING, MAXERRORSTRING);
}

void setErrorLoc(char errorloc[nErrorLoc]){
  strcopyN(ERROR_LOC, errorloc, nErrorLoc);
}

 

extern "C" {
#include <R_ext/Rdynload.h>

static R_NativePrimitiveArgType RelaxRFoption_t[] = { LGLSXP };
static const R_CMethodDef cMethods[]  = {
  {"RelaxUnknownRFoption", (DL_FUNC) &RelaxUnknownRFoption, 1, RelaxRFoption_t},
  {NULL, NULL, 0}
};


#define CALLDEF_DO(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  CALLDEF_DO(CholPosDef, 1),
  CALLDEF_DO(SolvePosDef, 3),
  CALLDEF_DO(struve, 4),
  CALLDEF_DO(I0ML0, 1),
  //  CALLDEF_DO(),
  {NULL, NULL, 0}
};


#define EXTDEF_DO(name, n)  {#name, (DL_FUNC) &name, n}
static const R_ExternalMethodDef extMethods[] = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  EXTDEF_DO(RFoptions, -1), 
  {NULL, NULL, 0} 
};


 

#define CALLABLE(FCTN)  R_RegisterCCallable("RandomFieldsUtils", #FCTN, (DL_FUNC)  FCTN)
void R_init_RandomFieldsUtils(DllInfo  *dll) {
  CALLABLE(solve_DELETE);
  CALLABLE(solve_NULL);
  CALLABLE(solvePosDef);
  CALLABLE(invertMatrix);
  
  CALLABLE(sqrtPosDef); 
  CALLABLE(sqrtRHS);

  CALLABLE(I0mL0);

  CALLABLE(getErrorString);
  CALLABLE(setErrorLoc);
  CALLABLE(getUtilsParam);
  CALLABLE(attachRFoptions);
  CALLABLE(detachRFoptions);
  CALLABLE(relaxUnknownRFoption);

  R_registerRoutines(dll, cMethods, callMethods, NULL, // .Fortran
		     extMethods);
  R_useDynamicSymbols(dll, FALSE);
}

void R_unload_RandomFieldsUtils(DllInfo VARIABLE_IS_NOT_USED *info) {
  /* Release resources. */
}

} // end extern C

/*

R
library("Matrix", lib="~/TMP/matrix")


 */
