#ifndef __BRENC_H__
#define __BRENC_H__

struct State;

#ifndef __INTPAIR_H__
#include "IntPair.h"
#endif

#ifndef __VECTOR_H__
#include "vector.h"
#endif

// Need the next one for NTAB.
#ifndef __BRENNER_H__
#include "brenner.h"
#endif

struct AtomPairInfo
{
  /* lcheck is 0, 1, or 2.  1 means both ends are either hydrogen or carbon.
     2 means both ends are either silicon or germanium.  0 means some other
     combination.  We set this up in the main caguts routine. */
  short lcheck;
  /* cor is the vector from the first atom to the second. */
  vector cor;
  /* rpp is the pairwise repulsion between them. */
  vector rpp;
  /* rcor is the length of cor. */
  Float rcor;
  /* ww is the value of f super c sub i j defined in equation 18 on page 16 for
     rcor. */
  Float ww;
  /* dww is the value of partial (f super c sub i j (r)) per partial r
     evaluated at rcor. */
  Float dww;
  /* exx1 is the attraction between the pair, V super A sub i j in equation 5
     on page 7. */
  Float exx1;
  /* dexx1 is the partial of exx1 with respect to rcor. */
  Float dexx1;
};

void mtable(Float drtable[4][4][NTAB],
            /* ddtab is the amount we multiply a radius by to get an index into
               tabfc and tabdfc.  Note that it used to be an amount we *divide*
               by, but I changed this  5 Aug 2000 because multiplication is
               faster.  Tim Freeman */
            Float ddtab[4][4],
            Float rtable[4][4][NTAB],
            Float atable[4][4][NTAB], Float datable[4][4][NTAB], 
            Float tabfc[4][4][NTAB], Float tabdfc[4][4][NTAB]);

void init_xh(FILE *inter2D);
void init_in2(FILE *inter2D);

#ifndef FUNGIMOL
void find_atom_energies(BrennerMainInfo *info);
#endif

Float RADIC(int KI, int KJ, Float XNT1, Float XNT2, Float CONJUG,
	     Float *drdl, Float *drdm, Float *drdn);
Float TOR( Float XNT1, Float XNT2, Float CONJUG,
	    Float *DATORI, Float *DATORJ, Float *DATORC);
Float BCUINT(int KJ, Float XX1, Float XX2, Float *ANSY1, Float *ANSY2);

#ifndef FUNGIMOL
int write_coord_file(FILE *fp, BrennerMainInfo *info);
#endif

void init_c ();
void deallocate_caguts ();
Float caguts(struct State *info);

struct AtomPairInfoState;
Float pibond(struct State *info, struct AtomPairInfoState *apis);
int choose_lj (struct State *s);
Float ljPure (struct State *s);
#ifndef FUNGIMOL
Float ljguts(BrennerMainInfo *state);
#endif
/* ljparam_init_kt should be called at initialization time to record the
   lennard jones parameters for each atom. kt1 is the ktype for the atom, which
   is kt[natom] where natom is the atomic number.  sigt and epst are the third
   and fourth columns of input.d. */
void ljparam_init_kt(int kt1, Float sigt, Float epst);
void ljparam(int ktmax, Float RLL);

void deallocate_rangetree ();
void make_neighbor_lists_tree(const struct State *s,
                              struct NeighborState *ljns,
                              Float rslj[NTYPES+1][NTYPES+1],
                              Float xmms[NTYPES+1][NTYPES+1]);

#ifndef FUNGIMOL
void vscale(BrennerMainInfo *info);
#endif
#endif /* BRENC_H */
