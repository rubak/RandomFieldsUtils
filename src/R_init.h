


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

#ifndef rfutils_init_H
#define rfutils_init_H 1


void solve_DELETE(solve_storage **S); 
void solve_NULL(solve_storage* x);

int solvePosDef_(double *M, int size, bool posdef, 
		double *rhs, int rhs_cols,
		double *logdet, 
		 solve_storage *PT, solve_param *sp,  int PL
		);

int invertMatrix(double *M, int size);

double I0mL0(double x);

void getErrorString(char errorstring[MAXERRORSTRING]);
void setErrorLoc(char errorloc[nErrorLoc]);


#endif

