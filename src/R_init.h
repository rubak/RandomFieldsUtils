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

