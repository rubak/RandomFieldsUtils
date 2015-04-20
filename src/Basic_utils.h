#ifndef basic_rfutils_h
#define basic_rfutils_h 1

#ifndef showfree
#define showfree !true 
#define DOPRINT true
// // 1
// #define RANDOMFIELDS_DEBUGGING 1
// // 1
#endif

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
