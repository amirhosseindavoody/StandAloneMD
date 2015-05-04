/*
 This code was released into the public domain by Tim Freeman on 15 Aug 2000.
*/

/* This is the C-to-C version of the api. It exists so I can replace it with a
   C-to-C++ version when I do the Fungimol integration. 
   Tim Freeman 15 Aug 2000.
*/

#ifndef __CCAPI_H__
#define __CCAPI_H__

#ifndef __API_H__
#include "api.h"
#endif

/* Need brenner.h for BrennerMainInfo. */
#ifndef __BRENNER_H__
#include "brenner.h"
#endif

/* A State is a state of the atoms in the system.  For the C to C interface,
   it's a BrennerMainInfo. */
struct State {
  BrennerMainInfo bmi;
};

static inline struct State *toState (BrennerMainInfo *bmi) {
  return (struct State *) bmi;
}

static inline const struct State *toStateconst (const BrennerMainInfo *bmi) {
  return (const struct State *) bmi;
}

static inline struct State *fromState (BrennerMainInfo *s) {
  return (struct State *) s;
}

/* Returns the first Atom. */
/* I'm having to go through and find all the loops from 0 to
   the number of atoms and make sure they start at firstIndex instead of 0,
   since stomping on atom 0 is not legitimate for the Fungimol plugin.  To
   avoid ever having to do this again in future code reorganizations, I'm
   making C-only code also use this function.  It's all inline so it should be
   just as fast. Tim Freeman 30 Oct 2000.
*/
static inline int firstIndex (const struct State *s) {
  return 0;
}

static inline int pastLastIndex (const struct State * s) {
  return s->bmi.num_atms;
}

static inline dvector pos (int ai, const struct State * s) {
  return s->bmi.atm_num[ai].coord;
}

static inline Double posx (int ai, const struct State * s) {
  return s->bmi.atm_num[ai].coord.x;
}

static inline Double posy (int ai, const struct State * s) {
  return s->bmi.atm_num[ai].coord.y;
}

static inline Double posz (int ai, const struct State * s) {
  return s->bmi.atm_num[ai].coord.z;
}

/* We need the previous coordinate so we know when to recompute the
   neighbors. */
static inline Double prevx (int ai, const struct State *s) {
  return s->bmi.atm_num[ai].prev_coord.x;
}

static inline Double prevy (int ai, const struct State *s) {
  return s->bmi.atm_num[ai].prev_coord.y;
}

static inline Double prevz (int ai, const struct State *s) {
  return s->bmi.atm_num[ai].prev_coord.z;
}

static inline void setPrevToCurrent (int ai, struct State *s) {
  BrenAtom *atm_num = s->bmi.atm_num;
  atm_num[ai].prev_coord.x = atm_num[ai].coord.x;
  atm_num[ai].prev_coord.y = atm_num[ai].coord.y;
  atm_num[ai].prev_coord.z = atm_num[ai].coord.z;
}

static inline Float rll (const struct State *s) {
  return s->bmi.RLL;
}

static inline void transForcex (int i, int j, struct State * s, Float f) {
  s->bmi.atm_num[i].force.x += f;
  s->bmi.atm_num[j].force.x -= f;
}

static inline void transForcey (int i, int j, struct State * s, Float f) {
  s->bmi.atm_num[i].force.y += f;
  s->bmi.atm_num[j].force.y -= f;
}

static inline void transForcez (int i, int j, struct State * s, Float f) {
  s->bmi.atm_num[i].force.z += f;

  s->bmi.atm_num[j].force.z -= f;
}

static inline int getKtype (int ai, const struct State * s) {
  return s->bmi.atm_num[ai].ktype;
}

static inline int getNoa (int ktype, const struct State *s) {
  return s->bmi.noa[ktype];
}

#ifndef INFINITE_CUBE
static inline const Float *getCube (const struct State * s) {
  return s->bmi.cube;
}
#endif

static inline int squeeze (const struct State *s) {
  return s->bmi.squeeze;
}

static inline int movable (int ai, const struct State *s) {
  return s->bmi.atm_num[ai].movable;
}

static inline struct NeighborState *ljNeighbors (struct State * s) {
  return s->bmi.ljNeighborState;
}

static inline const struct NeighborState *ljNeighborsconst (const struct State * s) {
  return s->bmi.ljNeighborState;
}

static inline struct NeighborState *caNeighbors (struct State * s) {
  return s->bmi.caNeighborState;
}

static inline const struct NeighborState *caNeighborsconst
(const struct State * s)
{
  return s->bmi.caNeighborState;
}
#endif
