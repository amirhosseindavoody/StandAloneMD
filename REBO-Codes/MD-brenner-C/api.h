/*
 This code was released into the public domain by Tim Freeman on 11 Aug 2000.
*/

#ifndef __API_H__
#define __API_H__

/* 

This is an API that allows C code to access atoms and neighbors of atoms.  The
goal is to allow one version of the source code for the C port of Brenner's
molecular dynamics code to either become a Fungimol plugin or to be part of a
purely C program that does batch manipulation of text files without any
significant loss of efficiency.  We get the required flexibility by having
multiple implementations of this api, and we get efficiency by having lots of
inline functions and being willing to recompile the library to use it in
multiple contexts rather than just relink it.

For starters, we just have the version of the API that supports the batch
program.

*/

/* Need vector.h because some of the API functions declared below use
   vectors. */
#ifndef __VECTOR_H__
#include "vector.h"
#endif

#ifndef __BOOL_H__
#include "bool.h"
#endif

struct NeighborState;
struct State;
static inline int firstIndex (const struct State *s);
static inline int pastLastIndex (const struct State *s);
static inline dvector pos (int ai, const struct State *s);
/* gcc doesn't optimize inline functions that return structures very well, so
   we have to provide inline functions that get the individual fields too. */
static inline Double posx (int ai, const struct State *s);
static inline Double posy (int ai, const struct State *s);
static inline Double posz (int ai, const struct State *s);
/* We need the previous coordinate so we know when to recompute the
   neighbors. */
static inline Double prevx (int ai, const struct State *s);
static inline Double prevy (int ai, const struct State *s);
static inline Double prevz (int ai, const struct State *s);
static inline void setPrevToCurrent (int ai, struct State *s);
/* rll governs how often we recompute neighbors.  We recompute neighbors if
   anybody moved rll/2 from their position last time we computed neighbors.
   The bigger rll is, the more non-neighbos will be included in each neighbor
   list.  rll is used to suitably inflate the neighbor cutoffs.  */
static inline Float rll (const struct State *s);
// The next three add the specified amount to the force for the first index
// and decrement for the second index.
static inline void transForcex (int i, int j, struct State *s, Float f);
static inline void transForcey (int i, int j, struct State *s, Float f);
static inline void transForcez (int i, int j, struct State *s, Float f);
static inline int getKtype (int ai, const struct State *s);
static inline int getNoa (int ktype, const struct State *s);
static inline int squeeze (const struct State *s);
static inline int movable (int ai, const struct State *s);
#ifndef INFINITE_CUBE
static inline const Float *getCube (const struct State *s);
#endif
static inline struct NeighborState *ljNeighbors (struct State *s);
static inline const struct NeighborState *ljNeighborsconst (const struct State *s);
static inline struct NeighborState *caNeighbors (struct State *s);
static inline const struct NeighborState *caNeighborsconst (const struct State *s);

void setAllPrevToCurrent (struct State *s);
#endif
