/*
 This code was spliced out of ljguts.c, which was released into the public
 domain by Peter McCluskey on 6/2/2000.  I'm leaving it in the public domain.
 Tim Freman 28 Oct 2000.
*/

/* For NTYPES. */
#ifndef __KT_H__
#include "kt.h"
#endif

/* This is the maximum cutoff for lennard-jones neighboring.  That is the
   maximum distance for which the force field has any effect, plus the slop so
   we don't have to reneighbor every time. */
extern Float rslj[NTYPES+1][NTYPES+1];
