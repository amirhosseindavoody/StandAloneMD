/* brenner.h - Copyright (c) 1998 Zyvex LLC.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    or its derived works must display the following acknowledgement:
 * 	This product includes software developed by Zyvex LLC.
 * 
 * This software is provided "as is" and any express or implied warranties,
 * including, but not limited to, the implied warranties of merchantability
 * or fitness for any particular purpose are disclaimed. In no event shall
 * Zyvex LLC be liable for any direct, indirect, incidental, special,
 * exemplary, or consequential damages (including, but not limited to,
 * procurement of substitute goods or services; loss of use, data, or
 * profits; or business interruption) however caused and on any theory of
 * liability, whether in contract, strict liability, or tort (including
 * negligence or otherwise) arising in any way out of the use of this
 * software, even if advised of the possibility of such damage.
 */

#include "kt.h"

#ifndef __BRENNER_H__
#include "brenner.h"
#endif

int ktmax;
int kt[MAX_ATOMNO+1];
int kt2[NTYPES+1];

void static_initkt () {
  int i;
  for(i = 0; i <= NTYPES; ++i) {
    kt2[i] = 0;
  }
  for (i = 0; i <= MAX_ATOMNO; i++) {
    kt[i] = 0;
  }
  kt[HYDROGEN] = 2;
  kt[CARBON] = 1;
  /* The next two used to write to kt out of bounds, back before I made it have
     the size MAX_ATOMNO+1 instead of NTYPES+1.  Tim Freeman  5 Oct 2000. */
  kt[SILICON] = 3;
  kt[GERMANIUM] = 4;
  kt2[2] = HYDROGEN;
  kt2[1] = CARBON;
  kt2[3] = SILICON;
  kt2[4] = GERMANIUM;
  ktmax = 4;
}
