


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
  void calcja_(int*, int*, int*, int*, int*, int*, int*);
  void amuxmat_(int*, int*, int*, double*, double*, double*, int*, int*);
  //  void transpose_(int *, int *, double *, int * int *, double*, int*, int*);
  //  void spamback_();
  //  void spamforward();
}


#endif
