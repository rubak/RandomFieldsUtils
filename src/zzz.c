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

#define RFU_NEED_OBSOLETE 1
#define NO_SSE2 1

#include "Basic_utils.h"
#include "win_linux_aux.h"
#include "RandomFieldsUtils.h"
#include "Utils.h"
#include "zzz_RandomFieldsUtils.h"

#define none 0

#if defined(__clang__)
//# pragma clang diagnostic ignored "-Wcast-function-type"
#endif

#ifdef __GNUC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
//// GCC diagnostic ignored "-Wcast-function-type"
#endif

static R_NativePrimitiveArgType 
    int_arg[] = { INTSXP },
    host_arg[] = { STRSXP, INTSXP};
  //  static R_NativeArgStyle argin[] = {R_ARG_IN},
  //    argout[] = {R_ARG_OUT},
  //   hostarg[] = {R_ARG_OUT, R_ARG_OUT};

#define CDEF(name, n, type) {#name, (DL_FUNC) & name, n, type}
static const R_CMethodDef cMethods[]  = {
  CDEF(sleepMilli,  1, int_arg),
  CDEF(sleepMicro, 1, int_arg),
  CDEF(pid, 1, int_arg),
  CDEF(hostname, 2, host_arg),
  CDEF(setCPUs, 1, int_arg),
  CDEF(recompilationNeeded, 1, int_arg),
  CDEF(loadoptions, 0, none),
  CDEF(detachoptions, 0, none),
  {NULL, NULL, 0, NULL}
};


#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  CALLDEF(SIMDmessages, 1),
  CALLDEF(DebugCall, 0),
  CALLDEF(Chol, 1),
  CALLDEF(debuggingLevel, 0),
  CALLDEF(scalarR, 3),
  CALLDEF(SolvePosDefR, 3),
  CALLDEF(struve, 4),
  CALLDEF(besselk_simd, 2),
  CALLDEF(I0ML0, 1),
  CALLDEF(gaussr, 2),
  CALLDEF(WMr, 4),
  CALLDEF(logWMr, 4),
  CALLDEF(sortX, 4),
  CALLDEF(orderX, 4), 
  CALLDEF(DivByRow, 2),
  CALLDEF(colMaxs, 1),
  CALLDEF(quadratic, 2),
  CALLDEF(dotXV, 2),
  CALLDEF(rowMeansX, 2),
  CALLDEF(rowProd, 1),
  CALLDEF(dbinorm, 2),
  CALLDEF(chol2mv, 2),
  CALLDEF(tcholRHS, 2),
  CALLDEF(crossprodX, 3),
  CALLDEF(getPackagesToBeInstalled, 1),
  CALLDEF(isGPUavailable,0),
  CALLDEF(isNEONavailable,0),
  CALLDEF(isX86_64,0),
  CALLDEF(gpu_info,1),
  CALLDEF(instruction_set, 3),
  //  CALLDEF(),
  {NULL, NULL, 0}
};


 
#define EXTDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_ExternalMethodDef extMethods[] = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  EXTDEF(RFoptions, -1), 
  {NULL, NULL, 0} 
};



#define CALLABLE(FCTN)  R_RegisterCCallable("RandomFieldsUtils", #FCTN, (DL_FUNC)  FCTN)

void R_init_RandomFieldsUtils(DllInfo  *dll) {
  CALLABLE(del_utilsoption);
  CALLABLE(get_utilsoption);
  CALLABLE(push_utilsoption);
  CALLABLE(params_utilsoption);

  CALLABLE(solve_DELETE);
  CALLABLE(solve_NULL);
  
  CALLABLE(SolvePosDef);
  CALLABLE(SolvePosDefSp);
  CALLABLE(SqrtPosDefFree);
  CALLABLE(xCinvXdet);
  CALLABLE(xCinvYdet);
  CALLABLE(DetPosDefsp);
  CALLABLE(InvertMatrix);
  CALLABLE(cholesky);
  CALLABLE(DetPosDef);
  CALLABLE(Is_positive_definite);

  

 
  CALLABLE(sqrtRHS); 
   CALLABLE(chol2inv);

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
  
  CALLABLE(attachRFUoptions);
  CALLABLE(detachRFUoptions);
  //  CALLABLE(linkRFUoptions);
  CALLABLE(RFUoptions);
  CALLABLE(attachSetNGet);
  CALLABLE(getoptionsRFU);
  CALLABLE(setoptionsRFU);

  
  // OBSOLETE_RFU
  CALLABLE(getUtilsParam);
  CALLABLE(attachRFoptions);
  CALLABLE(detachRFoptions);
  CALLABLE(relaxUnknownRFoption); // obsolete
  CALLABLE(getErrorString); // obsolete
  CALLABLE(setErrorLoc); // obsolete
  CALLABLE(ToIntI);
  CALLABLE(solvePosDef);
  CALLABLE(solvePosDefSp);
  CALLABLE(sqrtPosDefFree);
  CALLABLE(XCinvXdet);
  CALLABLE(XCinvYdet);
  CALLABLE(detPosDefsp);
  CALLABLE(detPosDef);
  CALLABLE(invertMatrix);
  CALLABLE(chol);
  CALLABLE(is_positive_definite);



  CALLABLE(ordering);
  CALLABLE(orderingL);
  CALLABLE(orderingInt);
  CALLABLE(orderingLong);
  CALLABLE(sorting);
  CALLABLE(sortingL);
  CALLABLE(sortingInt);
  CALLABLE(sortingLong);
  CALLABLE(scalarX);
  //  CALLABLE(scalarInt);

  CALLABLE(pid);
  CALLABLE(parallel);
  CALLABLE(sleepMicro); // problem?
 
  R_registerRoutines(dll, cMethods, callMethods, NULL, // .Fortran
		     extMethods);
  R_useDynamicSymbols(dll, FALSE); //
}


#ifdef SCHLATHERS_MACHINE
#ifdef __GNUC__ 
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
//// GCC diagnostic push
//// GCC diagnostic ignored "-Wcast-function-type"
#endif
#endif
void R_unload_RandomFieldsUtils(DllInfo *info) { }
#ifdef __GNUC__ 
//// GCC diagnostic pop
#endif
