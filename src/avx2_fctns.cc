
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of system specific auxiliary functions

 Copyright (C) 2021 -- 2021 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.
g
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include "Basic_utils_local.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include "RandomFieldsUtils.h"
#include "kleinkram.h"
#include "options.h"
#include "Utils.h"
#include "xport_import.h"
#include "extern.h"


#if defined AVX2

ASSERT_SIMD(avx2_fctns, avx2);

#define algn_general(X)  ((1U + (uintptr_t) (((uintptr_t) X - 1U) / BytesPerBlock)) * BytesPerBlock)

#if defined SSE41 || defined AVX2
int static inline *algnInt(int *X) {
  assert(algn_general(X)>=(uintptr_t)X); return (int *) algn_general(X);
}
#endif


void colMaxsIint256(int *M, Long r, Long c, int *ans) {
  if (r < 32
#if defined AVX2
      || !avx2Avail
#elif defined  SSE41
      || !sse41Avail
#endif      
       ) {
    for (int i=0; i<c; i++) {
      int *m = M + r * i,
	dummy = m[0];    
      for (int j=1; j<r; j++) dummy = MAX(dummy, m[j]);
      ans[i] = dummy;
    }
    return;
  }
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(static)
#endif    
  for (int i=0; i<c; i++) {
     int dummy,
      *m = M + r * i;
#if defined SSE41 || defined AVX2    
     int *start = algnInt(m),
       *end = m + r;
    uintptr_t End = (uintptr_t) (end - integers);
    if ((uintptr_t) start < End) {
      BlockType *m0 = (BlockType0*) start,
	Dummy = LOAD((BlockType0*) m0);
      for (m0++ ; (uintptr_t) m0 < End; m0++) {
	Dummy = MAXINTEGER(Dummy, LOAD(m0));
      }
      int *d = (int *) &Dummy;
      dummy = d[0];
      dummy = MAX(dummy, d[1]);
      dummy = MAX(dummy, d[2]);
      dummy = MAX(dummy, d[3]);
#if defined AVX2
      dummy = MAX(dummy, d[4]);
      dummy = MAX(dummy, d[5]);
      dummy = MAX(dummy, d[6]);
      dummy = MAX(dummy, d[7]);
#endif // AVX2
      for ( ; m<start; m++) dummy = MAX(dummy, *m);
      m = (int *) m0;
      for ( ; m<end; m++) dummy = MAX(dummy, *m);
    } else {
      dummy = m[0];    
      for (int j=1; j<r; j++) dummy = MAX(dummy, m[j]);
    }
#else // not SSE4
    dummy = m[0];    
    for (int j=1; j<r; j++) dummy = MAX(dummy, m[j]);
#endif    
    ans[i] = dummy;
  }
}


#else

void colMaxsIint(int *M, Long r, Long c, int *ans);
void colMaxsIint256(int *M, Long r, Long c, int *ans) {colMaxsIint(M, r, c, ans); }

SIMD_MISS(avx2_fctns, avx2);

#endif

