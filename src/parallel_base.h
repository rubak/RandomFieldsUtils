/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2021 -- 2021 Martin Schlather

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

#ifndef parallel_omp_base_H
#define parallel_omp_base_H 1

// #define NEVER_OMP 1
// #define NEVER_AVX 1
// #define NEVER_SSE 1

#if defined __WIN32__ || defined (__APPLE__) || defined(__sun)
  #if defined TIME_AVAILABLE
    #undef TIME_AVAILABLE
  #endif
#else
#define TIME_AVAILABLE 1
#endif


#if defined _OPENMP && ! defined NO_OMP && ! defined NEVER_OMP
  #if defined SCHLATHERS_MACHINE
    #define DO_PARALLEL 1
  #else
    #define DO_PARALLEL 1
  #endif
#elif defined DO_PARALLEL
  #undef DO_PARALLEL
#endif


#if defined NEVER_SSE
  #ifndef NO_SSE2
     #define NO_SSE2 1
  #endif
#elif defined NEVER_AVX
  #ifndef NO_AVX
     #define NO_AVX 1
  #endif
#endif


#if defined NO_SSE2 && ! defined NO_SSE3
#define NO_SSE3 1
#endif
#if defined NO_SSE3 && ! defined NO_SSSE3
#define NO_SSSE3 1
#endif
#if defined NO_SSSE3 && ! defined NO_SSE41
#define NO_SSE41 1
#endif
#if defined NO_SSE41 && ! defined NO_AVX
#define NO_AVX 1
#endif
#if defined NO_AVX && ! defined NO_AVX2
#define NO_AVX2 1
#endif
#if defined NO_AVX2 && ! defined NO_AVX512
#define NO_AVX512 1
#endif



#if defined AVX512
#undef AVX512
#endif

#if defined __AVX2__ && ! defined NO_AVX2
#define AVX2 1
#elif defined AVX2
#undef AVX2
#endif

#if defined __AVX__  && ! defined NO_AVX
#define AVX 1
#elif defined AVX
#undef AVX
#endif

#if defined __SSE41__ && ! defined NO_SSE41
#define SSE41 1
#elif defined SSE41
#undef SSE41
#endif

#if defined  __SSSE3__  && ! defined NO_SSSE3
#define SSSE3 1
#elif defined SSSE3
#undef SSSE3
#endif

#if defined  __SSE3__ && ! defined NO_SSE3
#define SSE3 1
#elif defined SSE3
#undef SSE3
#endif

#if defined  __SSE2__ && ! defined NO_SSE2
#define SSE2 1 
#elif defined SSE2
#undef SSE2
#endif



#endif
