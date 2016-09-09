
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of system specific auxiliary functions

 Copyright (C) 2001 -- 2015 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <math.h>
#include <unistd.h>
#include <Rinternals.h>
#include "win_linux_aux.h"
#include "General_utils.h"
#include "RandomFieldsUtils.h"

SEXP getChar() {
  ERR("does not work");
#ifdef WIN32
  ERR("input limitations on windows");
#endif
#define maxGetChar 255
  //typedef char intchar[sizeof(int) / sizeof(char)];
  //typedef char intchar[sizeof(int) / sizeof(char)];
  SEXP str;
  int //g,
    i = -1;
  char // c, *t = NULL,
    *s = NULL
    ;
  s = (char*) MALLOC(sizeof(char) * maxGetChar);
  // initscr();
  //  fflush(stdin);
  //  nocbreak();
  while (++i < maxGetChar) {
    // g = getchar();       
    //s[i] = ((intchar*) &g)[0][0];
    // g = scanf("%c\n", s); s[1] = '\0'; break;
    //t = fgets(s, 2, stdin); break;
    //    s[i] = getch();
    //fflush(stdin);
    if (false) {
      s[i+1] = '\0';
      // printf("%d i=%d  '%c' '%c' '%c' '%c' '%c'\n", g, i,  s[i],
      // ((intchar*) &g)[0][0],
      //  ((intchar*) &g)[0][1],
      // ((intchar*) &g)[0][2],
      // ((intchar*) &g)[0][3]
      // );
    }
    if (s[i] == '\n') {
      s[i] = '\0';
      break;
    }
  }
  //endwin();
//printf(">%s<\n", s);
  PROTECT(str=allocVector(STRSXP, 1));
  SET_STRING_ELT(str, 0, mkChar(s));  
  UNPROTECT(1);
  FREE(s);
  return str;
}
