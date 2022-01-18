
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of system specific auxiliary functions

 Copyright (C) 2001 -- 2021 Martin Schlather

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


#include "intrinsics.h"

#if defined MSDOS_WINDOWS
#define VC_EXTRALEAN
#include <windows.h>
#endif


// achtung! windows.h zusammen mit <Rmath.h oder R.graphics>
// gibt warnung, da ERROR mehrfach definiert !
// deshalb auch in auxiliary.h nicht basic.h einbinden // obsolette ?!!
#include <unistd.h>
#include <Rinternals.h>
#include "win_linux_aux.h"


 
void sleepMilli(int *milli) {
#if defined MSDOS_WINDOWS
  Sleep((long) *milli); // along
#else 
  usleep((useconds_t) (1000 * (unsigned long) *milli));// along
#endif
}

void sleepMicro(int *micro) {
#if defined MSDOS_WINDOWS
  Sleep((long) ((*micro + 500) / 1000));// along
#else
  usleep((useconds_t) *micro);
#endif
}

void pid(int *i)  {
#if defined MSDOS_WINDOWS
  *i = _getpid();
#else
  *i = getpid(); 
#endif
}

int parentpid=0;
bool
  parallel() {
    int mypid;
    pid(&mypid);
    return mypid != parentpid;
  }


void hostname(char **h, int *i){
#if defined MSDOS_WINDOWS
  *h[0]=0;
#else
  gethostname(*h, *i);
#endif
}  


uint32_t cpuid_info(int Blatt, int Register) {
#if defined MINGWCPUID
   uint32_t s[4];						     
   __cpuid(Blatt, s[0], s[1], s[2], s[3]);			    
   return s[Register];
#elif defined WINCPUID
  uint32_t s[4];							
  __cpuid((int *)s, (int) Blatt);
  return s[Register];
#elif defined LINUXCPUID
  uint32_t s[4];							
  asm volatile							
    ("cpuid": "=a"(s[0]), "=b"(s[1]),"=c"(s[2]),			
     "=d"(s[3]):"a"(Blatt),"c"(0));					
  return s[Register]; 
#else
   return 0;
#endif
}
