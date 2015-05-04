/* xmass.h - Copyright (c) 1998 Zyvex LLC.
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

#ifndef __XMASS_H__
#define __XMASS_H__

/* Need kt.h for NTYPES. */
#ifndef __KT_H__
#include "kt.h"
#endif

/* BOLZ is only used to compute ECONV. */
#define BOLZ 1.380662
/* EPSI0 is only used to compute ECONV. */
#define EPSI0 11604.5 /* EPSI as initialized in setin.f */
/* EPSI is only used in bren2.c when dealing with temperatures. */
#define EPSI 11605.0
/* AVO is only used to compute ECONV. */
#define AVO 6.02205
/* To convert the mass of an atom (as given in the xmass array) into real
   units, multiply by ECONV. */
#define ECONV ((1/(AVO*BOLZ*EPSI0))*1e7)

/* If an atom has atomic number n, then xmass [kt [n]] is the mass of one of
   them.  This is given a default value in alloc_bren in bren2.c, and it may
   be given a non-default value when reading input.d */
extern Float xmass[NTYPES+1];
/* emass is like xmass, except each element is multiplied by ECONV. */
extern Float emass[NTYPES+1];

/* You'd better call this before using xmass or emass. */
void static_initxmass ();
#endif
