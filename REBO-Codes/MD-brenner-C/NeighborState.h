/*
 This code was released into the public domain by Tim Freeman on 8/23/2000.
*/

/* This is a C structure for representing which atoms are neighbors of
   which. */
#ifndef __NEIGHBORSTATE_H__
#define __NEIGHBORSTATE_H__

struct NeighborState;

// FIXME I figured out a good conventions for ordering argument lists -- put
// the larger data structures earlier in the argument list.  By this
// convention, most of the argument lists here are wrong.  Fix it someday.
// It would be important to minimize elapsed time to minimize merge conflicts.
// Tim Freeman 13 Sep 2000.
int numNeighbors (int ai, const struct NeighborState *s);

const int *neighbors (int ai, const struct NeighborState *s);

// Make room for the given number of atoms. 
void allocate_atom (struct NeighborState *ljns, int max_atoms);

struct NeighborState *newNeighborState ();

void deallocateNeighborState (struct NeighborState *ljns);

static inline void allocate_neighbor (struct NeighborState *ljns,
                                      int lj_neighbor_end);

#endif

