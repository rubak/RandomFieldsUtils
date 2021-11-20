/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

 Copyright (C) 2019 -- 2021 Martin Schlather

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

#ifndef rfutils_intrinsics_H
#define rfutils_intrinsics_H 1

#ifdef PRINTF
#error "intrinsics.h not very first"
#endif

#include<inttypes.h> // uintptr_t
#include "def.h"
#include "parallel_simd.h"



#if defined WIN32 || defined _WIN32 || defined __WIN32__
//#warning loading intrin.h the first time
#include <intrin.h>							
#endif

#if defined AVX  || defined SSE2 //|| defined AVX2 || defined 
#include <immintrin.h>
#endif

#if __GNUC__ > 4 ||							\
  (__GNUC__ == 4 && (__GNUC_MINOR__ > 9 ||				\
		     (__GNUC_MINOR__ == 9 &&  __GNUC_PATCHLEVEL__ >= 1)))
//#define OpenMP4 1
#endif


union uni32{
  float f4[1];
  uint32_t u32[1];
  uint16_t u16[2];
  uint8_t u8[4];
};

union uni64{
  uint64_t u64;
  double d8;
  float f4[2];
  uint32_t u32[2];
  uint16_t u16[4];
  uint8_t u8[8];
};

union uni128{
#if defined SSE2
  __m128i vi; // hamming
  __m128d d;
  __m128 f;
  __m128d d128[1];
#endif
  uint64_t u64[2];// hamming
  uint32_t u32[4];
  uint8_t u8[16]; // hamming, shuffle
  //  __m64 m64[2];
  double halfd[2], d8[2];
  float halff[4], f4[4];
};


union uni256 {
#if defined AVX2
  __m256i vi; // hamming
  __m256d d;
  __m128d d128[2], halfd[2];
  __m256 f;
  __m128 halff[2];
#endif
  uint64_t u64[4];// hamming
  uint32_t u32[8];
  uint8_t u8[32]; // hamming, shuffle
  //  __m64 m64[4];
  double d8[4];
  float f4[8];
};


#define BitsPerByte 8L


#if defined AVX512
#define SSEBITS 512L
#define SSEMODE 30L
#define BIT_SHUFFLE _mm512_bitshuffle_epi64_mask
#define MASK0ADDFLOAT(A,M,B) A = _mm256_maskz_add_ps(A, M, A, B)


#elif defined AVX
#define SSEBITS 256L
#define SSEMODE 20L
#define BlockType0 __m256i 
#define BlockType __m256i ALIGNED
#define BlockUnitType0 uni256
#define Double __m256d
#define MAXDOUBLE _mm256_max_pd
#define LOAD _mm256_load_si256
// #define EXPDOUBLE mm256_exp_pd // only on intel compiler
#define ADDDOUBLE  _mm256_add_pd
#define SUBDOUBLE  _mm256_sub_pd
#define MULTDOUBLE _mm256_mul_pd 
#define LOADuDOUBLE _mm256_loadu_pd
#define LOADDOUBLE _mm256_load_pd
#define STOREuDOUBLE _mm256_storeu_pd
#define ZERODOUBLE _mm256_setzero_pd

#define MULTFLOAT  _mm256_mul_ps 
#define ADDFLOAT  _mm256_add_ps 
#define SUBFLOAT   _mm256_sub_ps 
#define ZEROFLOAT _mm256_setzero_ps
#define BLENDFLOAT  _mm256_blend_ps
#define DUPLICATEFLOAT  _mm256_moveldup_ps
#define MASK0ADDDOUBLE(A,M,B)  _mm256_maskz_add_pd(A, M, A, B)
#define BLENDDOUBLE  _mm256_blend_pd
#define DUPLICATEDOUBLE  _mm256_movedup_pd


#if defined AVX2
  #define MAXINTEGER _mm256_max_epi32

#define AND  _mm256_and_si256
#define OR  _mm256_or_si256
#define XOR  _mm256_xor_si256
#define ANY(A) (!_mm256_testz_si256(A, A))
#define SHR32  _mm256_srli_epi32 // see also _mm256512_rol_epi64,
#define SHL32  _mm256_slli_epi32
#define SHR16  _mm256_srli_epi16
#define SHR64  _mm256_srli_epi64
#define SHL64  _mm256_slli_epi64
#define SHUFFLE8(A,B,C) A=  _mm256_shuffle_epi8(B,C)

#define SET8  _mm256_set1_epi8
#define SETREV8(B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15)\
  _mm256_setr_epi8(B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,\
		   B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15)
#define SET16  _mm256_set1_epi16
#define SET32(A, B) A = _mm256_set1_epi32(B)
#define SET64  _mm256_set1_epi64x // oder _m256d _mm256_set1_pd (double a)
#define ZERO   _mm256_setzero_si256
#define LOAD  _mm256_load_si256
#define LOADU  _mm256_loadu_si256 // _mm256_lddqu_si256
#define STORE_DOUBLE _mm256_store_pd
#define EXTRACT16  _mm256_extract_epi16

#define ADD8 _mm256_add_epi8
#define ADD32  _mm256_add_epi32
#define MADD16  _mm256_madd_epi16 
#define ADD64  _mm256_add_epi64
#define SAD8   _mm256_sad_epu8
#define MULT32  _mm256_mullo_epi32

#else
  #define MAXINTEGER _mm_max_epi32
#endif


#elif defined SSE2
#define SSEBITS 128L
#define SSEMODE 10L
#define BlockType0 __m128i
#define BlockType __m128i ALIGNED
#define BlockUnitType0 uni128
#define Double __m128d
#define MAXDOUBLE _mm_max_pd
#define MAXINTEGER _mm_max_epi32
#define LOAD _mm_load_si128
#define EXPDOUBLE _mm_exp_pd // only on intel compiler
#define ADDDOUBLE  _mm_add_pd
#define SUBDOUBLE  _mm_sub_pd
#define MULTDOUBLE _mm_mul_pd 
#define LOADuDOUBLE _mm_loadu_pd
#define LOADDOUBLE _mm_load_pd
#define STOREuDOUBLE _mm_storeu_pd
#define ZERODOUBLE _mm_setzero_pd

#define MULTFLOAT  _mm_mul_ps 
#define ADDFLOAT  _mm_add_ps 
#define SUBFLOAT   _mm_sub_ps 
#define ZEROFLOAT _mm_setzero_ps
#define BLENDFLOAT  _mm_blend_ps
#define DUPLICATEFLOAT  _mm_moveldup_ps

#define AND  _mm_and_si128
#define OR  _mm_or_si128
#define XOR  _mm_xor_si128
bool any128(__m128i A);
#define ANY(A) any128(A)
#define SHR32  _mm_srli_epi32 // see also _mm512_rol_epi64,
#define SHL32  _mm_slli_epi32
#define SHR16  _mm_srli_epi16
#define SHR64  _mm_srli_epi64
#define SHL64  _mm_slli_epi64

#define SET8  _mm_set1_epi8
#define SETREV8 _mm_setr_epi8
#define SET16  _mm_set1_epi16
#define SET32(A, B)  A = _mm_set1_epi32(B)
#define SET64  _mm_broadcastq_epi64 
#define ZERO   _mm_setzero_si128
#define LOAD  _mm_load_si128
#define LOADU  _mm_loadu_si128
#define STORE_DOUBLE _mm_store_pd
#define EXTRACT16  _mm_extract_epi16

#define ADD8  _mm_add_epi8
#define ADD32  _mm_add_epi32
#define ADD64  _mm_add_epi64
#define MADD16  _mm_madd_epi16 
#define SAD8   _mm_sad_epu8
#define INT2FLOAT  _mm_cvtepi32_ps
#define INT2DOUBLE _mm_cvtpi32_pd // very expensive

#define BLENDDOUBLE  _mm_blend_pd
#define DUPLICATEDOUBLE  _mm_movedup_pd
//#define MOVEMASK _mm_movemask_ps
//#define BLEND _mm_blend_pd //see also _mm512_mask_inserti64x4_mm_insert_epi64


#if defined SSSE3 // within SSE2
#define SHUFFLE8(A, B,C)  A = _mm_shuffle_epi8(B,C)
#endif


#elif defined MMX || defined PlainInteger64 // particularly Bit23
#define SSEBITS 64L
#define SSEMODE 0L
#define BlockType0 uint64_t
#define BlockType BlockType0
#define BlockUnitType0 uni64

#define AND(B,C)  (B) & (C)
#define  OR(B,C)  (B) | (C)
#define XOR(B,C)  (B) xor (C) 
#define SHR64(B,C)  (B) >> (C)
#define SHL64(B,C)  (B) << (C)
#define SET32(A, B) { Uint *LoCaL = (Uint*) &(A); LoCaL[0] = B; LoCaL[1] = B; }
#define ADD64(B,C)  (B) + (C)
#define ZERO()  0L
#if defined MMX
  #define ADD8(B,C)  (BlockType0) _mm_add_pi8((__m64) B, (__m64) C)
#endif
#define ANY

  
#else  
#define SSEBITS 32L
#define SSEMODE 0L
#define BlockType0 uint32_t
#define BlockType BlockType0
#define BlockUnitType0 uni32
  #if defined PlainInteger32
    #define AND(B,C)  (B) & (C)
    #define  OR(B,C)  (B) | (C)
    #define XOR(B,C)  (B) xor (C) 
    #define SHR32(B,C)  (B) >> (C)
    #define SHL32(B,C)  (B) << (C)
    #define SET32(A, B) A = B
    #define ADD64(B,C) (B) + (C)
    #define ZERO()  0L
    #define ANY
  #else
    #if defined __GNUC__
      #warning No specification of any SIMD.
    #endif
  #endif
#endif

#define BytesPerBlock (SSEBITS / BitsPerByte)
#define ALIGNED __attribute__ ((aligned (BytesPerBlock)))
#define doubles (BytesPerBlock / 8)
#define integers (BytesPerBlock / 8)



#if defined AVX
#define SCALAR_DEFAULT SCALAR_NEARFMA
#else
#define SCALAR_DEFAULT SCALAR_BASE
#endif


#include "intrinsics_specific.h"

#endif

