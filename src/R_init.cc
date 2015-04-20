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
#include "utils.h"
#include "R_init.h"


extern "C" {
#include <R_ext/Rdynload.h>
  
#define CALLABLE(FCTN)  R_RegisterCCallable("RandomFieldsUtils", #FCTN, (DL_FUNC)  FCTN)

void R_init_RandomFieldsUtils(DllInfo VARIABLE_IS_NOT_USED *info) {
  CALLABLE(solve_DELETE);
  CALLABLE(solve_NULL);
  CALLABLE(solvePosDef_);
  CALLABLE(invertMatrix);

  CALLABLE(I0mL0);

  CALLABLE(getErrorString);
  CALLABLE(setErrorLoc);
}

void R_unload_RandomFieldsUtils(DllInfo VARIABLE_IS_NOT_USED *info) {
  /* Release resources. */
}

} // end extern C

/*

R
library("Matrix", lib="~/TMP/matrix")


 */
