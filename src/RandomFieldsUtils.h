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
