

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



#ifndef rfutils_error_H
#define rfutils_error_H 1


#define NOERROR 0                 
#define ERRORMEMORYALLOCATION 1 
#define ERRORFAILED 2      /* method didn't work for the specified parameters */
#define ERRORM 3           /* a single error message */
#define ERRORNOTPROGRAMMEDYET 4



#ifdef SCHLATHERS_MACHINE
#define ERRLINE PRINTF("(ERROR in %s, line %d)\n", __FILE__, __LINE__);
#else
#define ERRLINE 
#endif


#define LENMSG 250
#define MAXERRORSTRING 1000
#define nErrorLoc 1000
#define LENERRMSG 2000
extern char MSG[LENERRMSG], BUG_MSG[LENMSG], MSG2[LENERRMSG],
  ERRORSTRING[MAXERRORSTRING], ERRORSTRING_OK[MAXERRORSTRING],
  ERRORSTRING_WRONG[MAXERRORSTRING], ERROR_LOC[nErrorLoc];


#define RFERROR error
#define ERR(X) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X); RFERROR(MSG);}
#define ERR1(X, Y) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X); \
    sprintf(MSG2, MSG, Y);					 \
    RFERROR(MSG2);}
#define ERR2(X, Y, Z) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X);\
    sprintf(MSG2, MSG, Y, Z);					\
    RFERROR(MSG2);}
#define ERR3(X, Y, Z, A) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X);	\
    sprintf(MSG2, MSG, Y, Z, A);					\
    RFERROR(MSG2);}
#define ERR4(X, Y, Z, A, B) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X);	\
    sprintf(MSG2, MSG, Y, Z, A, B);					\
    RFERROR(MSG2);}
#define ERR5(X, Y, Z, A, B, C) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X); \
    sprintf(MSG2, MSG, Y, Z, A, B, C);					\
    RFERROR(MSG2);}
#define ERR6(X, Y, Z, A, B,C,D) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X); \
    sprintf(MSG2, MSG, Y, Z, A, B, C, D);				\
    RFERROR(MSG2);}
#define ERR7(X, Y, Z,A,B,C,D,E) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X); \
    sprintf(MSG2, MSG, Y, Z, A, B, C, D, E);				\
    RFERROR(MSG2);}
#define FERR(X) strcpy(ERRORSTRING, X); DEBUGINFOERR;
#define SERR(X) { FERR(X); return ERRORM;}
#define CERR(X) { FERR(X); err=ERRORM; continue;}
#define FERR1(X,Y) sprintf(ERRORSTRING, X, Y); DEBUGINFOERR;  
#define SERR1(X,Y) { FERR1(X, Y); return ERRORM;}
#define CERR1(X,Y) { FERR1(X, Y); err=ERRORM; continue; }
#define SERR2(X, Y, Z) { sprintf(ERRORSTRING, X, Y, Z);  DEBUGINFOERR; return ERRORM;}
#define CERR2(X, Y, Z) { sprintf(ERRORSTRING, X, Y, Z);  err=ERRORM; continue;}
#define SERR3(X, Y, Z, A) { sprintf(ERRORSTRING, X, Y, Z, A); DEBUGINFOERR; return ERRORM;}
#define CERR3(X, Y, Z, A) { sprintf(ERRORSTRING, X, Y, Z, A); err=ERRORM; continue;}
#define SERR4(X, Y, Z, A, B) { sprintf(ERRORSTRING, X, Y, Z, A, B); DEBUGINFOERR;  return ERRORM;}
#define SERR5(X, Y, Z, A, B, C) { sprintf(ERRORSTRING, X, Y, Z, A, B, C); DEBUGINFOERR; return ERRORM;}
#define SERR6(X, Y, Z, A, B, C, D) { sprintf(ERRORSTRING, X, Y, Z, A, B, C, D);  DEBUGINFOERR; return ERRORM;}
#define SERR7(X, Y, Z, A, B, C, D, E) { sprintf(ERRORSTRING, X, Y, Z, A, B, C, D, E); DEBUGINFOERR; return ERRORM;}
#define GERR(X) { strcpy(ERRORSTRING, X); err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR1(X,Y) { sprintf(ERRORSTRING, X, Y); err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR2(X,Y,Z) { sprintf(ERRORSTRING, X, Y, Z); err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR3(X,Y,Z,A) { sprintf(ERRORSTRING, X, Y, Z, A);  err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR4(X,Y,Z,A,B) { sprintf(ERRORSTRING, X, Y, Z, A, B);  err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR5(X,Y,Z,A,B,C) { sprintf(ERRORSTRING, X, Y, Z, A, B, C);  err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}

void getErrorString(char errorstring[MAXERRORSTRING]);
void setErrorLoc(char errorloc[nErrorLoc]);



#endif
