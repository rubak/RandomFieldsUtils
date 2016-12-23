


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
#include "General_utils.h"

extern "C" {
  void R_init_RandomFieldsUtils(DllInfo *info);
  void R_unload_RandomFieldsUtils(DllInfo *info);

  SEXP struve(SEXP X, SEXP Nu, SEXP Factor_Sign, SEXP Expscaled);
  SEXP I0ML0(SEXP X);
  SEXP gaussr(SEXP X, SEXP Derivative); 
  SEXP WMr(SEXP X, SEXP Nu, SEXP Derivative, SEXP Factor);
  SEXP logWMr(SEXP X, SEXP Nu1, SEXP Nu2, SEXP Factor);

  SEXP SolvePosDef(SEXP M, SEXP rhs, SEXP logdet);
  SEXP Chol(SEXP M);
  

  SEXP RFoptions(SEXP options);
  void RelaxUnknownRFoption(int *relax);

  SEXP attachRFoptionsUtils();
  SEXP detachRFoptionsUtils();

  SEXP sortX(SEXP Data, SEXP From, SEXP To, SEXP NAlast);
  SEXP orderX(SEXP Data, SEXP From, SEXP To, SEXP NAlast);

  void sleepMicro(int *micro);
  void sleepMilli(int *milli);
  void hostname(char **h, int *i);
  void pid(int *i);
  SEXP getChar();


  
  void Ordering(double *d, int *len, int *dim, int *pos);

}


#endif 
 
