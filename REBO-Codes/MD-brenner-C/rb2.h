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

#ifndef __RB2_H__
#define __RB2_H__
/* This used to be part of the BrennerMainInfo structure, but it's constant
   data and therefore logically static.  Also, the Fungimol plugin needs this
   data and it doesn't need the BrennerMainInfo structure.  Hence I'm making it
   static global.  Tim Freeman  4 Oct 2000. */
/* rb2[i][j] is D super max sub i j in equation 18 on page 16.  When we're
   computing f super c sub i j, and the radius is greater than D super max
   sub i j, then the f super c sub i j is 0. */ 
extern Float rb2[4][4];

#endif
