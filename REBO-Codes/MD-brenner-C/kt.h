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

#ifndef __KT_H__
#define __KT_H__

/* NTYPES: maximum # of different types of atoms allowed for Lennard-Jones 
 *        potential
 */
#define NTYPES 8

/* MAX_ATOMNO is the maximum atomic number, so we can size the kt table
   correctly.  I need a number that's big enough for both brennerc and
   everything Fungimol might ever do.  200 seems pretty safe.  Memory is
   cheap.
*/
#define MAX_ATOMNO 200

/* This is the maximum ktype, inclusive.  A ktype is a small internal index
   that uniquely identifies an atom type.  The kt2 array gives the mapping.
   We use elements 1 through ktmax inclusive of many arrays, including kt2,
   xmms, xm, rmaxlj, and rslj. */ 
extern int ktmax;
extern int kt[MAX_ATOMNO+1];
extern int kt2[NTYPES+1];

/* The caller had better call initkt before using any of the above
   variables.  There's another function called initkt, so call this one
   static_initkt. */
void static_initkt ();

#endif
