/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2017 Martin Schlather

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

#include "Basic_utils_local.h"
#include "RandomFieldsUtils.h"  // must be before anything else
#include "General_utils.h"
#include "own.h"
#include "zzz_RandomFieldsUtils.h"
#include "Utils.h"

extern const char * ownprefixlist[ownprefixN], **ownall[ownprefixN];
extern int ownallN[ownprefixN];

SEXP attachRandomFieldsUtils(SEXP show) {
  attachRFoptions(ownprefixlist, ownprefixN,
		  ownall, ownallN,
  		  setparameterUtils, NULL,
		  getparameterUtils,
		  delparameterUtils,
		  0, true);
  if (INTEGER(show)[0]) {
#ifdef DO_PARALLEL
   PRINTF("'RandomFieldsUtils' will use OMP\n");
#else
   PRINTF("'RandomFieldsUtils' will NOT use OMP\n");
#endif
  }
  return R_NilValue;
}

SEXP detachRandomFieldsUtils(){
  detachRFoptions(ownprefixlist, ownprefixN);
  freeGlobals();
  return R_NilValue;
}
