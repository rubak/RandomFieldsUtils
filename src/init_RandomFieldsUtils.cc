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
#include <R_ext/Rdynload.h>
#include "General_utils.h"
#include "init_RandomFieldsUtils.h"
#include "RandomFieldsUtils.h"
#include "own.h"


void getErrorString(errorstring_type errorstring){
  strcopyN(errorstring, ERRORSTRING, MAXERRORSTRING);
}

void setErrorLoc(errorloc_type errorloc){
  strcopyN(ERROR_LOC, errorloc, nErrorLoc);
}

 
SEXP attachRFoptionsUtils() {
  //  NList = 0;
  
  // printf("UTx %ld\n", (long) getUtilsParam);

  attachRFoptions(ownprefixlist, ownprefixN, ownall, ownallN,
  		  setparameterUtils, NULL, getparameterUtils);
  return R_NilValue;
}

SEXP detachRFoptionsUtils() {
  detachRFoptions(ownprefixlist, ownprefixN);
  return R_NilValue;
}

extern "C" {

#include <R_ext/Rdynload.h>

  static R_NativePrimitiveArgType Relax_t[] = { LGLSXP },
    int_arg[] = { INTSXP },
    host_arg[] = { STRSXP, INTSXP};
  static R_NativeArgStyle argin[] = {R_ARG_IN},
    argout[] = {R_ARG_OUT},
    hostarg[] = {R_ARG_OUT, R_ARG_OUT};
static const R_CMethodDef cMethods[]  = {
  {"RelaxUnknownRFoption", (DL_FUNC) &RelaxUnknownRFoption, 1, Relax_t,argin}, {"sleepMilli", (DL_FUNC) &sleepMilli, 1, int_arg, argin},
 {"sleepMicro", (DL_FUNC) &sleepMicro, 1, int_arg, argin},
 {"pid", (DL_FUNC) &pid, 1, int_arg, argout},
 {"hostname", (DL_FUNC) &hostname, 2, host_arg, hostarg},
  // {"attachRFoptionsUtils", (DL_FUNC) &attachRFoptionsUtils, 0, NULL, NULL},
  // {"detachRFoptionsUtils", (DL_FUNC) &detachRFoptionsUtils, 0, NULL, NULL},
  {NULL, NULL, 0, NULL, NULL}
};


#define CALLDEF_DO(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  CALLDEF_DO(Chol, 1),
  CALLDEF_DO(SolvePosDef, 3),
  CALLDEF_DO(struve, 4),
  CALLDEF_DO(I0ML0, 1),
  CALLDEF_DO(gaussr, 2),
  CALLDEF_DO(WMr, 4),
  CALLDEF_DO(logWMr, 4),
  CALLDEF_DO(attachRFoptionsUtils, 0),
  CALLDEF_DO(detachRFoptionsUtils, 0),
  CALLDEF_DO(sortX, 4),
  CALLDEF_DO(orderX, 4),
  CALLDEF_DO(getChar, 0), 
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
  CALLABLE(sqrtPosDefFree); 
  CALLABLE(sqrtRHS);

  CALLABLE(StruveH);
  CALLABLE(StruveL);
  CALLABLE(I0mL0);

  CALLABLE(WM);
  CALLABLE(DWM);
  CALLABLE(DDWM);
  CALLABLE(D3WM);
  CALLABLE(D4WM);
  CALLABLE(logWM);
  
  CALLABLE(Gauss);
  CALLABLE(DGauss);
  CALLABLE(DDGauss);
  CALLABLE(D3Gauss);
  CALLABLE(D4Gauss);
  CALLABLE(logGauss);
  

  CALLABLE(getErrorString);
  CALLABLE(setErrorLoc);
  CALLABLE(getUtilsParam);
  CALLABLE(attachRFoptions);
  CALLABLE(detachRFoptions);
  CALLABLE(relaxUnknownRFoption);

  CALLABLE(ordering);
  CALLABLE(orderingInt);
  CALLABLE(sorting);
  CALLABLE(sortingInt);

  R_registerRoutines(dll, cMethods, callMethods, NULL, // .Fortran
		     extMethods);
  R_useDynamicSymbols(dll, FALSE);
}

void R_unload_RandomFieldsUtils(DllInfo VARIABLE_IS_NOT_USED *info) {
  /* Release resources. */
}

}

/*

R
library("Matrix", lib="~/TMP/matrix")


 */
