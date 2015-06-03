


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


#ifndef RFutils_public_H
#define RFutils_public_H 1

#include <R_ext/Rdynload.h>
#include "utils.h"




extern "C" {
  void R_init_RandomFieldsUtils(DllInfo VARIABLE_IS_NOT_USED *info);
  void R_unload_RandomFieldsUtils(DllInfo VARIABLE_IS_NOT_USED *info);

  SEXP struve(SEXP X, SEXP Nu, SEXP Factor_Sign, SEXP Expscaled);
  SEXP I0ML0(SEXP X);
  SEXP solvePosDef(SEXP M, SEXP rhs, SEXP logdet);

  void R_init_RFutils(DllInfo *info);
}

#endif
