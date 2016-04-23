


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



#ifndef rfutils_matmult_H
#define rfutils_matmult_H 1


double XkCXtl(double *X, double *C, int nrow, int dim, int k, int l);
void XCXt(double *X, double *C, double *V, int nrow, int dim);
void AtA(double *a, int nrow, int ncol, double *A);
void xA(double *x, double*A, int nrow, int ncol, double *y);
void xA(double *x1, double *x2,  double*A, int nrow, int ncol, double *y1,
	double *y2);
void Ax(double *A, double*x, int nrow, int ncol, double *y);
void Ax(double *A, double*x1, double*x2, int nrow, int ncol, double *y1,
	double *y2);
double xUy(double *x, double *U, double *y, int dim);
double xUxz(double *x, double *U, int dim, double *z);
double x_UxPz(double *x, double *U, double *z, int dim);
double xUx(double *x, double *U, int dim);
void matmult(double *A, double *B, double *C, int l, int m, int n);
void matmulttransposed(double *A, double *B, double *C, int m, int l, int n);
void matmult_2ndtransp(double *A, double *B, double *C, int m, int l, int n);
void matmult_tt(double *A, double *B, double *C, int m, int l, int n);
double * matrixmult(double *m1, double *m2, int dim1, int dim2, int dim3);


#endif
