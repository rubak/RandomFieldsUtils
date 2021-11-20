/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 -- 2021 Martin Schlather

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
#include "Basic_utils.h"
#include "options.h"
#include "extern.h"
#include "zzz_RandomFieldsUtils.h"
#include "xport_import.h"
extern void strcopyN(char *dest, const char *src, int n);


#ifdef __cplusplus
extern "C" {
#endif

  typedef char errorstring_type[MAXERRORSTRING];
  bool RELAX_UNKNOWN_RFOPTION=false; // auf keinen Fall aendern!
  void relaxUnknownRFoption(bool relax){ RELAX_UNKNOWN_RFOPTION = relax; }
  
  void getErrorString(errorstring_type errorstring){
    STRCPY(errorstring, "error occurred in package RandomFieldsUtils");
  }
  
  void setErrorLoc(errorstring_type VARIABLE_IS_NOT_USED errorloc){
  }
  
#define MAXNLIST 7
#define PKGNAMELENGTH 20
  extern int NList;
  extern int noption_class_list,
    AllprefixN[MAXNLIST],
    *AllallN[MAXNLIST];
  extern  const char  *option_class_list[MAXNLIST],
    **Allprefix[MAXNLIST],
    ***Allall[MAXNLIST];
  extern  char pkgnames[MAXNLIST][PKGNAMELENGTH+1];
  extern setoptions_fctn setoption_fct_list[MAXNLIST];
  extern getoptions_fctn getoption_fct_list[MAXNLIST];
  extern finalsetoptions_fctn finaloption_fct_list[MAXNLIST];
  extern deleteoptions_fctn deloption_fct_list[MAXNLIST];
  extern bool installed [MAXNLIST];
  extern install_modes min_avx_needs[MAXNLIST],
    min_gpu_needs[MAXNLIST];
  extern Uint avx_infos [MAXNLIST];
  extern bool obsolete_package_in_use;

  setparameterfct setparam[MAXNLIST] = {NULL, NULL, NULL, NULL, NULL};
  getparameterfct getparam[MAXNLIST] = {NULL, NULL, NULL, NULL, NULL};
  finalsetparameterfct finalparam[MAXNLIST] = { NULL, NULL, NULL, NULL, NULL };
  deleteparameterfct delparam[MAXNLIST] = { NULL, NULL, NULL, NULL, NULL };
  
  
  void attachRFoptions(const char **PKGprefixlist, int N, 
		       const char ***PKGall, int *PKGallN,
		     setparameterfct set, finalsetparameterfct final,
		       getparameterfct get, deleteparameterfct del,
		       int pl_offset, bool basicopt) {
    char pkgname[] = "obsolete package";
    obsolete_package_in_use = true;
    RFU_GLOBAL_OPTIONS->solve.eigen2zero = 1e-10;
    RFU_GLOBAL_OPTIONS->basic.la_mode = 
      RFU_GLOBAL_OPTIONS->basic.la_usr = LA_INTERN;
     
    for (int ListNr=0; ListNr<NList; ListNr++) {    
      if (AllprefixN[ListNr] == N && 
	STRCMP(Allprefix[ListNr][0], PKGprefixlist[0]) == 0) {
	if (PL > 0) {
	  PRINTF("options starting with prefix '%.50s' have been already attached.",
		 PKGprefixlist[0]);
	}
	return;    
      }
    }
    if (basicopt) option_class_list[noption_class_list++] = PKGprefixlist[0];
    if (NList >= MAXNLIST) BUG;
    strcopyN(pkgnames[NList], pkgname, PKGNAMELENGTH);
    Allprefix[NList] = PKGprefixlist;
    AllprefixN[NList] = N;
    Allall[NList] = PKGall;
    AllallN[NList] = PKGallN;
    
    setoption_fct_list[NList] = NULL;
    finaloption_fct_list[NList] = NULL;
    getoption_fct_list[NList] = NULL;
    deloption_fct_list[NList] = NULL;
    
    setparam[NList] = set;
    finalparam[NList] = final;
    getparam[NList] = get;
    delparam[NList] = del;
    
    min_avx_needs[NList] = min_gpu_needs[NList] = Inone;
     
    NList++;
    PLoffset = pl_offset;
    PL = OPTIONS.basic.Cprintlevel = OPTIONS.basic.Rprintlevel + PLoffset;
    CORES = OPTIONS.basic.cores;
  }
  
  void detachRFoptions(const char **PKGprefixlist, int N) { detachRFUoptions(PKGprefixlist,  N); }

void getUtilsParam(utilsoption_type **global) { 
  *global = RFU_GLOBAL_OPTIONS; // OK!
}


  
#ifdef __cplusplus
}
#endif
