#ifndef basic_rfutils_h
#define basic_rfutils_h 1

#ifndef showfree
#define showfree !true 
#define DOPRINT true
// // 1
// #define RANDOMFIELDS_DEBUGGING 1
#endif


#define PL_IMPORTANT 1 
#define PL_SUBIMPORTANT 2
#define PL_RECURSIVE 3
#define PL_REC_DETAILS 4 
#define PL_STRUCTURE 5 // see also initNerror.ERROROUTOFMETHOD
#define PL_ERRORS  6 // only those that are caught internally

#define PL_FCTN_DETAILS 7  // R
#define PL_FCTN_SUBDETAILS 8

#define PL_COV_STRUCTURE 7 // C
#define PL_DIRECT_SEQU 8
#define PL_DETAILS 9
#define PL_SUBDETAILS 10


extern"C" {
  // Fortran Code by Reinhard Furrer
  void spamdnscsr_(int*, int*, double *, int*, double*, int*, int*, double*);
  void cholstepwise_(int*, int*, double*, int*, int*, int*, int*, int*,
		    int*, int*, int*, int*, int*, int*, double*, int*,
		    int*, int*, int*, int*);
  void backsolves_(int*, int*, int*, int*, int*, double*, int*, int*, int*,
		  int*, double*, double*);
  //  void transpose_(int *, int *, double *, int * int *, double*, int*, int*);
  //  void spamback_();
  //  void spamforward();
}


#endif
