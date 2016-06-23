
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 Martin Schlather

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

#include "RandomFieldsUtils.h"
#include "own.h"
#include "General_utils.h"
#include "kleinkram.h"
#include "own.h"




void setpDef(int VARIABLE_IS_NOT_USED  i, 
	     int VARIABLE_IS_NOT_USED  j, 
	     SEXP VARIABLE_IS_NOT_USED  el,
	     char VARIABLE_IS_NOT_USED  name[LEN_OPTIONNAME], 
	     bool VARIABLE_IS_NOT_USED  isList) {
  BUG;
}
void getpDef(SEXP VARIABLE_IS_NOT_USED  *sublist) {
  BUG;
}


bool RELAX_UNKNOWN_RFOPTION=false; // auf keinen Fall aendern!
void relaxUnknownRFoption(bool relax){ RELAX_UNKNOWN_RFOPTION = relax; }
void RelaxUnknownRFoption(int *relax){ relaxUnknownRFoption((bool) *relax); }




#define MAXNLIST 5
int NList = 1;
int AllprefixN[MAXNLIST] = {ownprefixN, 0, 0, 0, 0},
  *AllallN[MAXNLIST] = {ownallN, NULL, NULL, NULL, NULL};
const char **Allprefix[MAXNLIST] = {ownprefixlist, NULL, NULL, NULL, NULL},
  ***Allall[MAXNLIST] = { ownall, NULL, NULL, NULL, NULL};
setparameterfct setparam[MAXNLIST] = 
  {setparameterUtils, setpDef, setpDef, setpDef, setpDef};
getparameterfct getparam[MAXNLIST] = 
  {getparameterUtils, getpDef, getpDef, getpDef, getpDef};
finalsetparameterfct finalparam[MAXNLIST] = { NULL, NULL, NULL, NULL, NULL };

void setparameter(SEXP el, char *prefix, char *mainname, bool isList) {  
  int 
    j = NOMATCHING,
    i = NOMATCHING,
    ListNr = NOMATCHING;
  char name[LEN_OPTIONNAME];
  
  sprintf(name, "%s%s%s", prefix, strlen(prefix)==0 ? "" : ".", mainname);

  //  print("set param: %s.%s.%s\n",prefix, strlen(prefix)==0 ? "" : ".", mainname);
  //  print("relax=%d\n", RELAX_UNKNOWN_RFOPTION);

  if (mainname[0] >= 'A' && mainname[0] <= 'Z' && RELAX_UNKNOWN_RFOPTION) {
    if (PL > PL_IMPORTANT)
      PRINTF("'%s' is not considered as an RFoption, but will be passed to evaluate the model formula.\n", mainname);
    return;
  }

  if (strcmp(prefix, "")) {
    for (ListNr=0; ListNr<NList; ListNr++) {
      i = Match(prefix, Allprefix[ListNr], AllprefixN[ListNr]);
      // printf("Nr=%d i=%d \n", ListNr, i);
      if (i != NOMATCHING) break;
    }
    if (i == NOMATCHING) ERR1("option prefix name '%s' not found.", prefix); 
    if (i < 0 || strcmp(prefix, Allprefix[ListNr][i])) {
      for (int k=ListNr + 1; k < NList; k++) {
	int ii = Match(prefix, Allprefix[ListNr], AllprefixN[ListNr]);
	if (ii == NOMATCHING) continue;
	i = MULTIPLEMATCHING;
	if (ii >= 0 && strcmp(prefix, Allprefix[k][ii])==0) {
	  ListNr = k;
	  i = ii;
	  break;
	} // ii >0
      } // for k
      //printf("ii=%d %d %d\n", ii, NOMATCHING, i);
      if (i == MULTIPLEMATCHING) 
	ERR1("option prefix name '%s' is ambiguous.", prefix);   
    } // prefix == List


    // printf("ListNr = %d %s %d %s\n", ListNr, Allprefix[ListNr][i], i,
    //	   mainname);

    j = Match(mainname, Allall[ListNr][i], AllallN[ListNr][i]);
      
  } else { // (i==0), no prefix given
#define MinNameLength 3
    for (ListNr=0; ListNr<NList; ListNr++) {
      int trueprefixN = AllprefixN[ListNr];
      for (i=0; i < trueprefixN; i++) {	
	j = Match(mainname, Allall[ListNr][i], AllallN[ListNr][i]);
	if (j != NOMATCHING) break;     
      }
      if (j != NOMATCHING) break;
    }
    if (j==NOMATCHING) ERR1("Unknown option '%s'.", name);
    
    // printf("j=%d  %s\n", j,  j >=0 ? Allall[ListNr][i][j] : "multi");
    
    if (j < 0  || strcmp(mainname, Allall[ListNr][i][j])) {
      int starti = i + 1;
      for (int k = ListNr; k<NList; k++, starti=0) {
	int tpN = AllprefixN[k];
	for (int ii=starti; ii < tpN; ii++) {	
	  int jj = Match(mainname, Allall[k][ii], AllallN[k][ii]);
	  if (jj == NOMATCHING) continue;
	  
	  //  printf("listnr=%d %s jj=%d %s\n", ListNr, Allall[ListNr][i][j], 
	  //	 jj, jj < 0 ? "none" : Allall[k][ii][jj]);
	  j = MULTIPLEMATCHING;
	  if (jj >= 0 && strcmp(mainname, Allall[k][ii][jj])==0) {
	    ListNr = k;
	    i = ii;
	    j = jj;
	    break;
	  } // jj
	} // for ii
	if (j >= 0) break;
      } // for k
    } // if j < 0 || !=
  } // no prefix given

  if (j<0) ERR1("Multiple matching for '%s'.", name); 
 

  //  printf("%s %d %d %d \n", name, ListNr, i, j);
  
  setparam[ListNr](i, j, el, name, isList); 
}




SEXP getRFoptions() {
  SEXP list, names;
  
  int i, ListNr, itot = 0,
    k = 0;
  int  trueprefixN, totalN;

  for (totalN=ListNr=0; ListNr<NList; ListNr++) {
    trueprefixN = AllprefixN[ListNr];
    for (i=0; i<trueprefixN; i++) {
      //printf("ListNr=%d i=%d %s\n", ListNr, i, Allprefix[ListNr][i]);
      totalN += strcmp(Allprefix[ListNr][i], OBSOLETENAME) != 0;
    }
  }

  PROTECT(list = allocVector(VECSXP, totalN));
  PROTECT(names = allocVector(STRSXP, totalN));

  SEXP *sublist, *subnames;
  sublist = (SEXP *) MALLOC(sizeof(SEXP) * totalN);
  subnames = (SEXP *) MALLOC(sizeof(SEXP) * totalN);
  for (ListNr =0; ListNr<NList; ListNr++) {
    trueprefixN = AllprefixN[ListNr];    
    for (i=0; i<trueprefixN; i++, itot++) {   
      if (strcmp(Allprefix[ListNr][i], OBSOLETENAME) == 0) continue;
      SET_STRING_ELT(names, itot, mkChar(Allprefix[ListNr][i]));
      PROTECT(sublist[itot] = allocVector(VECSXP, AllallN[ListNr][i]));
      PROTECT(subnames[itot] = allocVector(STRSXP, AllallN[ListNr][i]));
      int endfor = AllallN[ListNr][i];
      for (k=0; k<endfor; k++) {
	SET_STRING_ELT(subnames[itot], k, mkChar(Allall[ListNr][i][k])); 
      }
    }
    getparam[ListNr](sublist + itot - trueprefixN);
  }

  for (i=0; i<totalN; i++) {
    setAttrib(sublist[i], R_NamesSymbol, subnames[i]);
    SET_VECTOR_ELT(list, i, sublist[i]);
  }
  setAttrib(list, R_NamesSymbol, names);
   
  assert(N == totalN);
  UNPROTECT(2 + 2 * totalN);
  FREE(sublist);
  FREE(subnames);

  return list;
}


void splitAndSet(SEXP el, char *name, bool isList) {
  int i, len;
  char prefix[LEN_OPTIONNAME / 2], mainname[LEN_OPTIONNAME / 2];   
  //  printf("splitandset\n");
  len = strlen(name);
  for (i=0; i < len && name[i]!='.'; i++);
  if (i==0) { ERR1("argument '%s' not valid\n", name); }
  if (i==len) {
    strcpy(prefix, "");
    strcopyN(mainname, name, LEN_OPTIONNAME / 2);
  } else {
    strcopyN(prefix, name, MIN(i + 1, LEN_OPTIONNAME / 2));
    strcopyN(mainname, name+i+1, MIN(strlen(name) - i, LEN_OPTIONNAME / 2) );
  }

  // 
  //printf("i=%d %d %s %s\n", i, len, prefix, mainname);
  setparameter(el, prefix, mainname, isList && GLOBAL.basic.asList);
  //   printf("ende\n");
}


SEXP RFoptions(SEXP options) {
  int i, j, lenlist, lensub;
  SEXP el, list, sublist, names, subnames;
  char *name, *pref;
  bool isList = false;
  /* 
     In case of strange values of a parameter, undelete
     the comment for PRINTF
  */

  
  // PRINTF("start %f\n", GLOBAL.gauss.exactness);
  options = CDR(options); /* skip 'name' */
  if (options == R_NilValue) {
    //PRINTF("before get %f\n", 1.);
    return getRFoptions(); 
  }

  name = (char*) (isNull(TAG(options)) ? "" : CHAR(PRINTNAME(TAG(options))));
  if ((isList = strcmp(name, "LIST")==0)) {   
    list = CAR(options);
    if (TYPEOF(list) != VECSXP)
      ERR1("'LIST' needs as argument the output of '%s'", RFOPTIONS);
    names = getAttrib(list, R_NamesSymbol);   
    lenlist = length(list);
    for (i=0; i<lenlist; i++) {
      int len;
      pref = (char*) CHAR(STRING_ELT(names, i));  

      //   print("%d %s warn.ambig=%d\n", i, pref, GLOBAL.warn.ambiguous);

      sublist = VECTOR_ELT(list, i);
      len = strlen(pref);
      for (j=0; j < len && pref[j]!='.'; j++);
      if (TYPEOF(sublist) == VECSXP && j==len) { // no "."
	// so, general parameters may not be lists,
	// others yes
	lensub = length(sublist);
	subnames = getAttrib(sublist, R_NamesSymbol); 
	for (j=0; j<lensub; j++) {
	  name = (char*) CHAR(STRING_ELT(subnames, j));

	  // print("    %d %s warn.ambig=%d\n", j, name, GLOBAL.warn.ambiguous);

	  //    print("%d %d %s : %f %f\n", i, j, name,
	  //	     GLOBAL.gauss.exactness, GLOBAL.TBM.linesimustep);	  
	  //
	  //	  print("  %d %d pref=%s name=%s\n", i, j, pref, name);
	  	  
	  setparameter(VECTOR_ELT(sublist, j), pref, name, 
		       isList & GLOBAL.basic.asList);
	}
      } else {   
	splitAndSet(sublist, pref, isList);
      }
    }
    //    print("end1 %f\n", GLOBAL.TBM.linesimufactor);
  } else {
    for(i = 0; options != R_NilValue; i++, options = CDR(options)) {

      //printf("set opt i=%d\n", i);

      el = CAR(options);
      name = (char*) (isNull(TAG(options)) ? "" :CHAR(PRINTNAME(TAG(options))));
      splitAndSet(el, name, isList);
    }
    //       print("end2 %f\n", GLOBAL.gauss.exactness);
  }

  for (i=0; i<NList; i++) if (finalparam[i] != NULL) finalparam[i]();

  GLOBAL.basic.asList = true;
  return(R_NilValue);
} 
 


void attachRFoptions(const char **prefixlist, int N, 
		  const char ***all, int *allN,
		  setparameterfct set, finalsetparameterfct final,
		  getparameterfct get) {
  for (int ListNr=0; ListNr<NList; ListNr++) {    
    if (AllprefixN[ListNr] == N && 
	strcmp(Allprefix[ListNr][0], prefixlist[0]) == 0) {
      if (PL > 0) 
	PRINTF("options starting with prefix '%s' have been already attached.",
		prefixlist[0]);
      return;    
    }
  }
  if (NList >= MAXNLIST) BUG;
  Allprefix[NList] = prefixlist;
  AllprefixN[NList] = N;
  Allall[NList] = all;
  AllallN[NList] = allN;
  setparam[NList] = set;
  finalparam[NList] = final;
  getparam[NList] = get;
  NList++;
}



void detachRFoptions(const char **prefixlist, int N) {
  int ListNr;
  for (ListNr=0; ListNr<NList; ListNr++) {    
    if (AllprefixN[ListNr] == N && 
	strcmp(Allprefix[ListNr][0], prefixlist[0]) == 0) break;
  }  
  if (ListNr >= NList) {
    ERR1("options starting with prefix '%s' have been already attached.",
	prefixlist[0]);
  }

  for (ListNr++; ListNr<NList; ListNr++) {    
    Allprefix[ListNr - 1] = Allprefix[ListNr];
    AllprefixN[ListNr - 1] =  AllprefixN[ListNr];
    Allall[ListNr - 1] = Allall[ListNr];
    AllallN[ListNr - 1] = AllallN[ListNr];
    setparam[ListNr - 1] = setparam[ListNr];
    finalparam[ListNr - 1] = finalparam[ListNr];
    getparam[ListNr - 1] = getparam[ListNr];
  }

  NList--;
}

void getUtilsParam(utilsparam **global) { 
  *global = &GLOBAL; 
}


