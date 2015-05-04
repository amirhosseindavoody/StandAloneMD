/* bren.c - Copyright (c) 1998 Zyvex LLC.
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

/* Reads in coordinates of atoms in Hyperchem, determines their
   energies, steps them in space to lower their energies, and updates
   Hyperchem coordinates.  Also maybe accept keyboard input for velocities
   of and forces on atoms.  Uses Brenner's potential: Phys. Rev. B (1990) v42 n15 p9458-70.
   Based on new version from http://www.engr.ncsu.edu/mat/CompMatSci/projects.html
   File HCnewpot1.ps is unpublished Brenner paper describing the new version. */

#ifndef __RLIST_H__
#include "rlist.h"
#endif

#ifndef __XMASS_H_
#include "xmass.h"
#endif

#ifndef __RMAX_H__
#include "rmax.h"
#endif

#ifndef __KT_H__
#include "kt.h"
#endif

#include "brenner.h"

#include "sili_germ.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "myassert.h"
#include <time.h>
#include <string.h>
#include "spgch.h"

#ifndef __BOOL_H__
#include "bool.h"
#endif

#ifndef __ANYAPI_H__
#include "anyapi.h"
#endif

#ifndef __CALC_DIST_H__
#include "calc_dist.h"
#endif

#ifndef __CCNEIGHBORSTATE_H__
#include "CCNeighborState.h"
#endif

#ifndef __RB2_H__
#include "rb2.h"
#endif

/* xalloc.h can define free as a macro, and free is called below. */
#ifndef __XALLOC_H__
#include "xalloc.h"
#endif

typedef void VOID;
#define RBOXES 64
/*#define BOXWIDTH 1.7 */
#define BOXWIDTH 4
#define BOXATMS 80
#define JITTER 1.e-3

/*Debye temperature converted to fs-1 X 2Pi*/
#define WD (2230*(2.08365e-5)*2*M_PI)

/* Fungimol can't afford the space for the boxes.  It doesn't use them anyway. */
#ifndef FUNGIMOL
/* A box has a list of atoms in it.  There are RBOXES^3 of them. */

struct BOX { int list[BOXATMS]; int numAtoms;};
typedef struct BOX SB;
static SB box[RBOXES][RBOXES][RBOXES];

/* The code below reads this from the third line of input.d.  (It used to read
   it from the second line, but that makes things incompatible with the sample
   input.d that came with the Fortran version, so I revised it to read from the
   third line.)  The comment at this place in the sample input.d from
   Brennermd's Execute directory says this is "=1 REBO (C,H,Si,Ge), =2
   tight-binding for C".  I disagree with Peter's comment (which I have
   deleted) saying that it's the periodic box length.  The dimensions of the
   periodic box are the fourth line of coord.d.  Tim Freeman 20 Jul 2000.*/

static int IPOT;

static int most_atoms_in_a_box;

static void initkt(BrennerMainInfo *info, int *kflag);
void append_file(FILE *fp, Float TTIME, Float DELTA, Float *cube);
static void initialize_atoms(BrennerMainInfo *info, int kflag);
static void apply_thermostat(BrennerMainInfo *);
static void third_order_corrector(BrennerMainInfo *);
static void third_order_predictor(BrennerMainInfo *);
static void make_neighbor_lists(BrennerMainInfo *);
static int look_in_this(SB *box1, BrenAtom *pAtm, BrennerMainInfo *state,
                        int kend);
static void find_boxes(BrennerMainInfo *);
#endif
/* Get coefficients for H and F and torsion. */
void PARAM(int ktmax, Float RLL)
{
  /*Parameters for hydrocarbons*/
  int I,IC,I2D,J,K,L,M,N;
  int i,j;
  Float XHH;
  FILE *inter2D, *inter3Dh, *inter3Dch;
  
  if ((inter2D = fopen(DATA_HOME "/inter2d_iv.d","r"))==NULL) {
    fprintf(stderr, "opening inter2d_iv failed\n");
    my_exit(-1);
  }
  fscanf(inter2D, "%d\n",&I2D);
  
  /*Read bicubic spline coefficients*/
  init_xh(inter2D);
  init_in2(inter2D);
  fclose(inter2D);
  for(i = 0; i < 4; ++i)
    {
      for(j = 0; j < 4; ++j)
        {
          Float r2 = rb2[i][j];
          RLIST[i+1][j+1] = (r2+RLL)*(r2+RLL);
#ifndef NDEBUG
          if(0) printf("rlist[%d][%d] %9.4f %9.4f %9.4f rmax %7.3f\n",
                       i,j,RLIST[i+1][j+1],r2,RLL, r2*r2);
#endif
          rmax[i][j] = r2*r2;
        }
    }
  ljparam(ktmax, RLL);
}

#ifndef FUNGIMOL

static int cmp_n(const void *v1, const void *v2)
{
  const IntPair *i1 = (const IntPair*)v1;
  const IntPair *i2 = (const IntPair*)v2;
  return i1->v1 == i2->v1 ? (i1->v2 - i2->v2) : (i1->v1 - i2->v1);
}
 
/* Compute neighbor lists if there are less than 200 atoms.  Straightforward
   n-squared algorithm, where we look at each atom and then look at all
   other atoms as possible neighbors.  The return value will be the number of
   neighbors that we found. */
static int
n2_neighbor_list(BrennerMainInfo *info)
{
  struct NeighborState *cans = info->caNeighborState;
  int num_atms = info->num_atms;
  BrenAtom *atm_num = info->atm_num;
  /* We will return k. */
  int k = 0;
  int i, j;
  allocate_atom (cans, num_atms);
  for(i = firstIndex (toState (info)); i < num_atms; ++i) {
    const Float *rlist_ki;
    int ki = atm_num[i].ktype;
    allocate_neighbor (cans, k+num_atms);
    cans->neighbor_start[i] = k;
    
    /*
     * cuts out all but C,H,Si, and Ge
     */
    if(ki >= 5) continue;
    rlist_ki = RLIST[ki];
    
    /* We're not making use of the fact that if atom i is a neighbor of
       j, then j is a neighbor of i.  This loop could probably start at j=i+1
       and run twice as fast. Tim Freeman 23 Jul 2000. */
    /* Added code to do that if N2_ASYM defined; the speed difference is
       too small to measure for typical molecules. pcm 2000-08-11 */
    /* And the new neighbor data structure has to be created in order, so the
       code Peter put in doesn't work any more and it isn't obvious how to fix
       it.  Thus I'm ripping out the code I badgered Peter into putting in
       here.  Sorry.  Tim Freeman 28 Aug 2000. */
    for(j = firstIndex (toState (info)); j < num_atms; ++j) {
      int kj;
      Float rlis, rsq;
      if(i == j) continue;
      kj = atm_num[j].ktype;
      /*
       * cuts out all but C,H,Si, and Ge
       */
      if(kj >= 5) continue;
      rlis = rlist_ki[kj];
      rsq = CALC_DIST_NO_RR(toState (info), i, j);
      if(rsq > rlis) continue;
      cans->neighbor_list[k] = j;
      ++k;
    }
  }
  cans->neighbor_start[i] = k;
  return k;
}

  static SB *boxes_used[MAX_ATOMS];
  static int num_boxes_used;

static void add_to_box(int x, int y, int z, int i, const BrenAtom *const atm_num)
{
    assert (0 <= z && z < RBOXES);
    if (box[x][y][z].numAtoms < BOXATMS) {
      box[x][y][z].list[box[x][y][z].numAtoms] = i;
      if(!box[x][y][z].numAtoms)
        boxes_used[num_boxes_used++] = &box[x][y][z];
      if(++box[x][y][z].numAtoms > most_atoms_in_a_box) {
        most_atoms_in_a_box = box[x][y][z].numAtoms;
        if(most_atoms_in_a_box > 48)
          printf("most_atoms_in_a_box %d\n", most_atoms_in_a_box);
      }
    } else {
      fprintf(stderr,
              "Too many atoms in box %d at %d, %d, %d\n",
              box[x][y][z].numAtoms, x, y, z);
      fprintf(stderr, "atom %d coords %.2f %.2f %.2f\n",
              i, atm_num[i].coord.x, atm_num[i].coord.y,
              atm_num[i].coord.z);
#ifdef INFINITE_CUBE
      fprintf(stderr,"This may indicate a need to increase RBOXES\n");
#endif
      my_exit(-1);
    }
}

/* Sort the atoms into boxes.  The boxes are in the file-scope array "box"
   declared above.
*/
static void find_boxes(BrennerMainInfo *info)
{
  int i, j, k;
  int x,y,z;
  int num_atms = info->num_atms;
  const BrenAtom *const atm_num = info->atm_num;
#ifndef INFINITE_CUBE
  const Float *cube = info->cube;
  const double cubex = cube[0]/BOXWIDTH;
  const double cubey = cube[1]/BOXWIDTH;
  const double cubez = cube[2]/BOXWIDTH;
  const int cubexmin = RBOXES/2 - (int)ceil(cubex/2);
  const int cubexmax = RBOXES/2 + (int)floor(cubex/2);
  const int cubeymin = RBOXES/2 - (int)ceil(cubey/2);
  const int cubeymax = RBOXES/2 + (int)floor(cubey/2);
  const int cubezmin = RBOXES/2 - (int)ceil(cubez/2);
  const int cubezmax = RBOXES/2 + (int)floor(cubez/2);
#endif
  
  /* initialize boxes to empty */
  for(i = 0; i < num_boxes_used; ++i)
    boxes_used[i]->numAtoms = 0;
  /* assign atoms to boxes */ 
  num_boxes_used = 0;
  for (i=firstIndex (toState (info)); i<num_atms; i++) {
    Double t;
    /* move the atoms +RBOXES/2 in every coordinate to put the atoms
       in the middle of the box group */
    t = atm_num[i].coord.x;
#ifndef INFINITE_CUBE
    t -= cube[0]*floor(t/cube[0] + 0.5);
#endif
    x = (int)floor(t/BOXWIDTH)+RBOXES/2;
    assert (0 <= x && x < RBOXES);
    t = atm_num[i].coord.y;
#ifndef INFINITE_CUBE
    t -= cube[1]*floor(t/cube[1] + 0.5);
#endif
    y = (int)floor(t/BOXWIDTH)+RBOXES/2;
    assert (0 <= y && y < RBOXES);
    t = atm_num[i].coord.z;
#ifndef INFINITE_CUBE
    t -= cube[2]*floor(t/cube[2] + 0.5);
#endif
    z = (int)floor(t/BOXWIDTH)+RBOXES/2;
    add_to_box(x, y, z, i, atm_num);
#ifndef INFINITE_CUBE
		/* add duplicate "wrapped around" entries */
    if (x == cubexmax)
      add_to_box(cubexmin, y, z, i, atm_num);
    else if (x == cubexmin)
      add_to_box(cubexmax, y, z, i, atm_num);
    if (y == cubeymax)
      add_to_box(x, cubeymin, y, i, atm_num);
    else if (y == cubeymin)
      add_to_box(x, cubeymax, y, i, atm_num);
    if (z == cubezmax)
      add_to_box(x, y, cubezmin, i, atm_num);
    else if (z == cubezmin)
      add_to_box(x, y, cubezmax, i, atm_num);
#endif
  }
}

/* Given an atom and a box, read box.list[] into the neighbor structure.
   The caller should have allocated enough room in the neighbor structure for
   this already.
   Return the new number of neighbors.
*/
static inline int
look_in_this(SB *box1, BrenAtom *pAtm, BrennerMainInfo *state, int kend)
{
	int i;
	int ki = pAtm->ktype;
	const Float *rlist_ki = RLIST[ki];
    const BrenAtom *const atm_num = state->atm_num;
	for (i = 0; i < box1->numAtoms; i++) {
		int index2 = box1->list[i];
		const BrenAtom *pAtm2 = &atm_num[index2];
		int kj = pAtm2->ktype;
		Float rlis, rsq;
		/*
		 * cuts out all but C,H,Si, and Ge
		 */
		if(kj >= 5) continue;
        /* FIXME We could put them in the neighbor list regardless of distance,
           and filter out the bogus ones when we scan the neighbor list later.
           This would make for longer neighbor lists, with the benefit of one
           less distance computation.  I don't know whether this is a net win.
           Peter says he'll replace this with a range tree anyway, so I'm not
           interested in trying this optimization soon.  Tim Freeman 28 Aug
           2000.
           But the code that would filter out the neighbors that are too far
           apart is already present in caguts.c, so the only cost of the
           experiment should be ifdef'ing out the code here.  Maybe worth
           trying anyway once I've transformed everything to use the API.  Tim
           Freeman 28 Aug 2000. */
		rlis = rlist_ki[kj];
		rsq = CALC_DIST_NO_RR(toState (state), pAtm->number, index2);
        /* If it's too far away, leave it out. */
		if(rsq > rlis) continue;
        /* The following if statement ensures that no atom is counted ase its
           own neighbor. */
		if (index2 != pAtm->number) {
          /* FIXME I hope the compiler figures out that
             state->caNeighborState->neighbor_list is constant and hoists it
             out of this inline code.  Seems dubious.  Might want to check the
             assembly language.
             Tim Freeman 28 Aug 2000 */
          state->caNeighborState->neighbor_list [kend] = index2;
          ++kend;
		}
	}
	return kend;
}

int
out_of_box(double x)
{
  return fabs(x) > RBOXES*BOXWIDTH/2;
}
/* Fill in neighbor_list to reflect who is a neighbor of who.  At the end, for
   each atom (suppose it's number is i), info->atom_num[i].num_neighbors will
   be the number of neighbors of atom i.  The neighbors will be represented as
   consecutive elements of neighbor_list.  The first one will be neighbor_list
   [info->atom_num[i].neighbor_start].  If j is a neighbor of i, then one
   element of neighbor_list will have v1 equal to i and v2 equal to j.  In the
   likely event that i is also a neighbor of j, then another element of
   neighbor_list will apparently have v1 equal to j and v2 equal to i. 
   Tim Freeman 24 Jul 2000. */
static void make_neighbor_lists(BrennerMainInfo *info)
{
	int i,j,k,n,x,y,z;
	int kend = 0;
	BrenAtom *pAtom;
	int num_atms = info->num_atms;
	BrenAtom *atm_num = info->atm_num;
    struct NeighborState *cans = info->caNeighborState;
#ifndef INFINITE_CUBE
	const Float *cube = info->cube;
#endif
	if(num_atms >= info->min_range_tree
#ifndef INFINITE_CUBE 
       || out_of_box(info->cube[0])
	   || out_of_box(info->cube[1])
       || out_of_box(info->cube[2])
#endif
	   )
	{
	  static Float xmms[NTYPES+1][NTYPES+1]; /* initialized by the compiler to
                                                all zeros */
	  make_neighbor_lists_tree(toState (info), info->caNeighborState,
                               RLIST, xmms);
	} else if(num_atms < 50) {
      n2_neighbor_list(info);
      printf("n2_neighbor_list\n");
    } else {
#ifndef INFINITE_CUBE
      const double cubex = info->cube[0];
      const double cubey = info->cube[1];
      const double cubez = info->cube[2];
#endif
      find_boxes(info);
      printf("bins\n");
      /*find the neighbor list for each atom by going through 
        the 27 surrouding boxes. */
      allocate_atom (cans, num_atms);
      for (n=firstIndex (toState (info)); n<num_atms; n++) {
        int xplus1, yplus1, zplus1;
        int xminus1, yminus1, zminus1;
        double t;
        allocate_neighbor (cans, kend + num_atms);
        cans->neighbor_start [n] = kend;
        pAtom = &atm_num[n];
        if(pAtom->ktype >= 5) continue;
        t = atm_num[n].coord.x;
#ifndef INFINITE_CUBE
        t -= cube[0]*floor(t/cube[0] + 0.5);
#endif
        x = (int)floor(t/BOXWIDTH)+RBOXES/2;
        /* Peter's code here used to care what the cube was, even in the infinite
           cube case.  I don't know the right thing to do in that case.  I think
           Peter's code was wrong in the infinite cube case because xminus1 might
           be negative if the x coordinate has become sufficiently negative, thus
           causing us to reference our array out of bounds.  Instead I'll use
           simpler code that I absolutely know is wrong, by just assuming the
           entire scene will fit within the grid of boxes.  This code will simply
           do the wrong thing if atoms in the scene escape from the region that
           we've made boxes for.  Tim Freeman 26 Aug 2000. */
#ifdef INFINITE_CUBE
        xminus1 = x - 1;
        xplus1 = x + 1;
#else
	xminus1 = x - 1;
	if (xminus1 < RBOXES/2 - cubex/(2*BOXWIDTH))
	  xminus1 = (int)floor(xminus1 + cubex/BOXWIDTH + 0.5);
	xplus1 = x + 1;
	if (xplus1 > RBOXES/2 + cubex/(2*BOXWIDTH))
	  xplus1 = (int)floor(xplus1 - cubex/BOXWIDTH + 0.5);
#endif
        t = atm_num[n].coord.y;
#ifndef INFINITE_CUBE
        t -= cube[1]*floor(t/cube[1] + 0.5);
#endif
        y = (int)floor(t/BOXWIDTH)+RBOXES/2;
#ifdef INFINITE_CUBE
        yminus1 = y - 1;
        yplus1 = y + 1;
#else
	yminus1 = y - 1;
	if (yminus1 < RBOXES/2 - cubey/(2*BOXWIDTH))
	  yminus1 = (int)floor(yminus1 + cubey/BOXWIDTH + 0.5);
	yplus1 = y + 1;
	if (yplus1 > RBOXES/2 + cubey/(2*BOXWIDTH))
	  yplus1 = (int)floor(yplus1 - cubey/BOXWIDTH + 0.5);
#endif
        t = atm_num[n].coord.z;
#ifndef INFINITE_CUBE
        t -= cube[2]*floor(t/cube[2] + 0.5);
#endif
        z = (int)floor(t/BOXWIDTH)+RBOXES/2;
#ifdef INFINITE_CUBE
        zminus1 = z - 1;
        zplus1 = z + 1;
#else
	zminus1 = z - 1;
	if (zminus1 < RBOXES/2 - cubez/(2*BOXWIDTH))
	  zminus1 = (int)floor(zminus1 + cubez/BOXWIDTH + 0.5);
	zplus1 = z + 1;
	if (zplus1 > RBOXES/2 + cubez/(2*BOXWIDTH))
	  zplus1 = (int)floor(zplus1 - cubez/BOXWIDTH + 0.5);
#endif
        kend = look_in_this(&box[xminus1][yminus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yminus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yminus1][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][y][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][y][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][y][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yplus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yplus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yplus1][zplus1],	pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yminus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yminus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yminus1][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][y][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][y][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][y][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yplus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yplus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yplus1][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yminus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yminus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yminus1][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][y][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][y][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][y][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yplus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yplus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yplus1][zplus1], pAtom,
                            info, kend);
      }
      cans->neighbor_start [n] = kend;
    }
}

static void third_order_predictor(BrennerMainInfo *info)
{
  int i;
  int num_atms = info->num_atms;
  BrenAtom *atm_num = info->atm_num;
  
  for (i=firstIndex (toState (info)); i<num_atms; i++) 
    if (atm_num[i].movable) {
      atm_num[i].coord =
        dplus(dplus(dplus(atm_num[i].coord,atm_num[i].velocity),
                    atm_num[i].accel),
              atm_num[i].dx3dt3);
      atm_num[i].velocity =
        plus(plus(atm_num[i].velocity,product(2,atm_num[i].accel)),
             product(3,atm_num[i].dx3dt3));
      atm_num[i].accel = plus(atm_num[i].accel,product(3,atm_num[i].dx3dt3));
    }						
}

static void third_order_corrector(BrennerMainInfo *info)
{
  int i;
  Float DE;
  vector RI;
  int num_atms = info->num_atms;
  BrenAtom *atm_num = info->atm_num;
  /* We don't really need to recompute DE each time.  However, doing it once
     per timestep is so small that the extra computation won't be
     measurable. */
  DE=-(info->timestep*info->timestep/2)/ECONV;
  
  for (i=firstIndex (toState (info)); i<num_atms; i++) 
    if (atm_num[i].movable) {
      RI = plus(atm_num[i].accel,product(DE/atm_num[i].mass,atm_num[i].force));
      atm_num[i].coord = dplus(atm_num[i].coord, product(-1.0/6,RI));
#ifndef NDEBUG
      if(!(atm_num[i].coord.x > -1.e9))
        {
          fprintf(stderr, "absurd coord %d %e xforce %e a %e mass %e\n",
                  i, RI.x, atm_num[i].force.x, atm_num[i].accel.x, atm_num[i].mass);
          my_exit(-1);
        }
#endif
      atm_num[i].velocity = plus(atm_num[i].velocity, product(-5.0/6,RI));
      /* Note that the old value of atm_num[i].accel was one of the summands in
         RI, so adding -RI to accel is a complicated way of ignoring the old
         value of accel.  This is good, since the force field really does tell
         us the new acceleration. */ 
      atm_num[i].accel = plus(atm_num[i].accel, product(-1.0,RI));
      atm_num[i].dx3dt3 = plus(atm_num[i].dx3dt3, product(-1/3.0,RI));
    }
}

/* FRICTION AND RANDOM FORCE */
/* I'm uncertain whether this works as well as the Fortran version; the
 * overall behavior looks suspiscious, but looking at individual changes
 * I don't see anything that can be explained by randomness. pcm 2000-08-30.
 */
static void
apply_gleq_thermostat(BrennerMainInfo *info)
{
  int i;
  int ii;
  int nta = 0;
  int nlr;
  const double PI2 = M_PI*2;
  const Float TR = info->temperature/EPSI/ECONV;
  const Float BET = WD*M_PI*ECONV/6.0/info->timestep;
  const Float GSIG = sqrt(2.0*TR*ECONV*BET);
  Float *gl;
  int *nlist = (int *)malloc(info->num_atms * sizeof(*nlist));

  for(i = firstIndex (toState (info)); i < info->num_atms; ++i) 
    if(info->atm_num[i].thermostated)
    {
      nlist[nta++] = i;
    }
  nlr = 3*nta/2;
  gl = (Float *)malloc(2*nlr*sizeof(gl[0]));

  for(i = 0; i < nlr; ++i)
  {
    double rr = rand()/(double)RAND_MAX;
    if(rr >= 1.0e-6)
    {
      float pre = sqrt(-2.0*log(rr));
      gl[i] = pre*GSIG*(Float)cos(PI2*rand()/RAND_MAX);
      gl[i+nlr] = pre*GSIG*(Float)sin(PI2*rand()/RAND_MAX);
    }
  }

  for(ii = 0; ii < nta; ++ii)
  {
    i = nlist[ii];
    if(info->atm_num[i].thermostated)
    {
      Float mass = xmass[info->atm_num[i].ktype];
      Float bm = BET*mass;
      Float sm = sqrt(mass);
      info->atm_num[i].force.x -= bm*info->atm_num[i].velocity.x + sm*gl[ii];
      info->atm_num[i].force.y -= bm*info->atm_num[i].velocity.y + sm*gl[ii+nta];
      info->atm_num[i].force.z -= bm*info->atm_num[i].velocity.z + sm*gl[ii+2*nta];
    }
  }
  xfree(gl);
  xfree(nlist);
}

/* USE EVANS-HOOVER SCHEME */

static void
apply_hoov_thermostat(BrennerMainInfo *info)
{
  /* this used for all atoms */
  int i;
  Float sc;
  Float ff = 0.0;
  Float df = 0.0;
  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    Float mass = xmass[info->atm_num[i].ktype];
    Float velx = info->atm_num[i].velocity.x;
    Float vely = info->atm_num[i].velocity.y;
    Float velz = info->atm_num[i].velocity.z;
    ff += info->atm_num[i].force.x*velx;
    df += velx*velx*mass;
    ff += info->atm_num[i].force.y*vely;
    df += vely*vely*mass;
    ff += info->atm_num[i].force.z*velz;
    df += velz*velz*mass;
  }
  if(df == 0)
  {
    fprintf(stderr, "Error in hoov thermostat - atom velocities are all zero?\n");
    my_exit(-1);
  }
  sc = ff/df;
  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    Float mass = xmass[info->atm_num[i].ktype];
    info->atm_num[i].force.x -= sc*info->atm_num[i].velocity.x*mass;
    info->atm_num[i].force.y -= sc*info->atm_num[i].velocity.y*mass;
    info->atm_num[i].force.z -= sc*info->atm_num[i].velocity.z*mass;
  }
}

static void apply_thermostat(BrennerMainInfo *info)	/*USE BERENDSEN SCHEME*/
{
	int i;
	Float XX, SC, SM;
	Float TR = info->temperature/EPSI/ECONV;
	Float BET = (WD*M_PI*ECONV/6/info->timestep);
	int num_atms = info->num_atms;
	BrenAtom *atm_num = info->atm_num;

	XX=0;
	for (i=firstIndex (toState (info)); i<num_atms; i++) 
		if (atm_num[i].thermostated)
			XX += dot(atm_num[i].velocity,atm_num[i].velocity)*atm_num[i].mass;

	if (XX<1e-7) 
	{
		fprintf(stderr, "T=0, Reset Thermostat to other than 1, fatal");
		exit(-1);
	}

	SC=BET*(TR*info->dnla/XX-1);

	for (i=firstIndex (toState (info)); i<num_atms; i++) 
		if (atm_num[i].thermostated) {
			SM=atm_num[i].mass*SC;
			atm_num[i].force = plus(atm_num[i].force,product(SM,atm_num[i].velocity));
			if(0)
			  printf("therm %2d %9.6f V %10.8f\n",
			       i, atm_num[i].force.x, atm_num[i].velocity.x);
		}
}

/* Despite the name, this finds forces too.  This is not static because it is
   also called from minimize.c. */
void
find_atom_energies(BrennerMainInfo *info)
{
  int i;
  int num_atms = info->num_atms;
  BrenAtom *atm_num = info->atm_num;
  static int cnt1;
  
  if(info->lchk == 1) {
    make_neighbor_lists(info);
#ifndef NDEBUG
    /* We can't validate ljNeighbors yet, since it's set up in ljguts. */
    validate (caNeighbors (toState(info)), 0, num_atms);
    sortNeighbors (info->caNeighborState, 0, num_atms);
    validate (caNeighbors (toState(info)), 0, num_atms);
    {
      int i;
      struct NeighborState *cans = caNeighbors (toState(info));
      for (i = firstIndex (toState(info)); i < num_atms; i++) {
	const int iNeighborCount = numNeighbors (i, cans);
	const int *const iNeighbors = neighbors (i, cans);
	int jn;
#if 0
	printf("atom %2d %2d nbors\n", i, iNeighborCount);
	for (jn = 0; jn < iNeighborCount; jn++) {
	  const int j = iNeighbors [jn];
	  printf("  %d\n", j);
	}
#endif
      }
    }    
#endif
  }
  info->system_energy = 0;
  for (i=firstIndex(toState(info)); i<num_atms; i++) {
    atm_num[i].force.x = 0;
    atm_num[i].force.y = 0;
    atm_num[i].force.z = 0;
  }
  
  info->system_energy += caguts(toState (info));
  info->system_energy += ljguts(info);
}

int
write_coord_file(FILE *fp, BrennerMainInfo *info)
{
  int i;
  const BrenAtom *atm_num = info->atm_num;
  fprintf(fp, "%s\n", info->head);
  fprintf(fp, "%6d\n", info->num_atms);
  fprintf(fp, "%20.11e %20.11e\n", info->starttime, info->timestep);
#ifdef INFINITE_CUBE
  fprintf(fp, "%20.11e %20.11e %20.11e\n", 1000.0, 1000.0, 1000.0);
#else
  fprintf(fp, "%20.11e %20.11e %20.11e\n", info->cube[0], info->cube[1],
	  info->cube[2]);
#endif
  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    fprintf(fp,"%5d %5d %20.11e %20.11e %20.11e %3d\n",
	    i+1, kt2[atm_num[i].ktype], atm_num[i].coord.x, 
	    atm_num[i].coord.y, atm_num[i].coord.z,
	    atm_num[i].movable ? atm_num[i].thermostated : 2);
  }

  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    fprintf(fp, "%5d %20.11e %20.11e %20.11e\n", i+1, atm_num[i].velocity.x, 
	    atm_num[i].velocity.y, atm_num[i].velocity.z);
  }

  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    fprintf(fp, "%5d %20.11e %20.11e %20.11e\n", i+1, atm_num[i].accel.x, 
	    atm_num[i].accel.y, atm_num[i].accel.z);
  }

  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    fprintf(fp, "%5d %20.11e %20.11e %20.11e\n", i+1, atm_num[i].dx3dt3.x, 
	    atm_num[i].dx3dt3.y, atm_num[i].dx3dt3.z);
  }
  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    fprintf(fp, "%5d %20.11e %20.11e %20.11e\n", i+1, 0.0, 0.0, 0.0);
  }
  return 1;
}

static void initialize_atoms(BrennerMainInfo *info, int kflag)
{
	int i;
	Float small;
	int num_atms = info->num_atms;
	BrenAtom *atm_num = info->atm_num;
    
    /* FIXME This value for small won't work in single-precision. */
	small=1e-12;
	for (i=firstIndex (toState (info)); i<num_atms; i++) 
      if (atm_num[i].movable==1){
        if(kflag != 4) {
          /* Add some noise, so that when we try to do the thermostat we
             never have zero temperature, since that leads to dividing by
             zero. */
          atm_num[i].velocity.x = (rand()/(Float)RAND_MAX)*2*JITTER-JITTER;
          atm_num[i].velocity.y = (rand()/(Float)RAND_MAX)*2*JITTER-JITTER;
          atm_num[i].velocity.z = (rand()/(Float)RAND_MAX)*2*JITTER-JITTER;
          atm_num[i].accel.x = (rand()/(Float)RAND_MAX)*2*small-small;
          atm_num[i].accel.y = (rand()/(Float)RAND_MAX)*2*small-small;
          atm_num[i].accel.z = (rand()/(Float)RAND_MAX)*2*small-small;
          
          atm_num[i].dx3dt3.x = (rand()/(Float)RAND_MAX)*2*small-small;
          atm_num[i].dx3dt3.y = (rand()/(Float)RAND_MAX)*2*small-small;
          atm_num[i].dx3dt3.z = (rand()/(Float)RAND_MAX)*2*small-small;
        }
        atm_num[i].prev_coord.x = atm_num[i].coord.x;
        atm_num[i].prev_coord.y = atm_num[i].coord.y;
        atm_num[i].prev_coord.z = atm_num[i].coord.z;
      }
}

#ifndef __SAFE_FGETS_H__
#include "safe_fgets.h"
#endif

static void
initkt (BrennerMainInfo *info, int *kflag)
{
  char buf[128];
  int i;
  int kuc, maxkb;
  FILE *fp = fopen("input.d", "r");
  if(!fp)
  {
    fprintf(stderr, "can't read input.d\n");
    my_exit(-1);
  }
  safe_fgets(buf, sizeof(buf), fp, "input.d line 1");
  if(sscanf(buf, "%d %d %d %d", &kuc, &maxkb, kflag, &info->nxmol) != 4)
  {
    fprintf(stderr, "parse error in input.d line 1\n");
    my_exit(-1);
  }
  /* Random # seed, neighbor list, temperature(K) */
  safe_fgets(buf, sizeof(buf), fp, "input.d line 2");
  /* I want to be able to feed the same files to the distributed Fortran
     brennermd program and to this program, so we need compatibility.  This
     code used to demand that ipot be on the same line as the rest.  Note that
     we never use IPOT.  Tim Freeman 24 Jul 2000*/
  {
    double pseed, rll, temp;
    if(sscanf(buf, "%lf %lf %lf", &pseed, &rll, &temp) != 3) {
      fprintf(stderr, "parse error in input.d line 2\n");
      my_exit(-1);
    }
    info->PSEED = pseed;
    info->RLL = rll;
    info->temperature = temp;
  }
  safe_fgets(buf, sizeof(buf), fp, "unused potential flag");
  if(sscanf(buf, "%d", &IPOT) != 1) {
    fprintf(stderr, "parse error in input.d line 3\n");
    my_exit(-1);
  }
  while(fgets(buf, sizeof(buf), fp)) {
      int natom;
      Float xma, epst, sigt;
      {
        double xmad, epstd, sigtd;
        if (sscanf(buf, "%d %lf %lf %lf", &natom, &xmad, &epstd , &sigtd) != 4) {
          fprintf (stderr, "Parse error in Lennard-Jones paramters for this "
                   "line:\n%s\n", buf);
          my_exit (-1);
        }
        xma = xmad;
        epst = epstd;
        sigt = sigtd;
      }
      if(natom < 0)
        continue;
      if(kt[natom] == 0) {
        if(ktmax > NTYPES) {
          fprintf(stderr, "Error - number of types (%d) exceeds limit\n",
                  ktmax);
          exit(-1);
        }
        kt[natom] = ktmax;
        kt2[ktmax] = natom;
      }
      xmass[kt[natom]] = xma;
      ljparam_init_kt(kt[natom], sigt, epst);
      ++ktmax;
  }
  fclose(fp);
}

static void
remove_com_velocity(BrennerMainInfo *info)
{
  Float xmt = 0.0;
  int i;
  vector com;
/*
****IMPORTANT********************************
*                                           *
* IF NO RIGID ATOMS, SUBTRACT COM VELOCITY  *
*                                           *
* REMOVE THIS CODE FOR MOLECULAR COLLISION  *
*                                           *
*********************************************
*/
  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    if(!info->atm_num[i].movable)
      return;
    xmt += xmass[info->atm_num[i].ktype];
  }
  com.x = com.y = com.z = 0.0;
  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    com.x += info->atm_num[i].velocity.x*xmass[info->atm_num[i].ktype];
    com.y += info->atm_num[i].velocity.y*xmass[info->atm_num[i].ktype];
    com.z += info->atm_num[i].velocity.z*xmass[info->atm_num[i].ktype];
  }
  com.x /= xmt;
  com.y /= xmt;
  com.z /= xmt;
  for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
  {
    info->atm_num[i].velocity.x -= com.x;
    info->atm_num[i].velocity.y -= com.y;
    info->atm_num[i].velocity.z -= com.z;
  }
}

#ifdef RUNGEKUTTA
/* Fourth order Runge-Kutta is a really bad choice for numerical integration
   here, because (a) It's fourth-order, and we don't have a bound on the fifth
   derivative of our result, so we don't have any error bound whatsoever and
   (b) it's a single step method, so we evaluate our force field many more
   times than needed.  However, it's what I have in Fungimol right now, so with
   a sufficiently small step size I hope to get good results that I'll be able
   to repeat after making this a Fungimol plugin.
   I don't want to use Nordsieck with Fungimol because the version of the
   algorithm stated in Brenner's code appears to require that we're computing
   position given acceleration, and that isn't the case for my DesignAtoms, and
   I think DesignAtoms are really cool so I don't want to break them.  I
   haven't found the original paper by Nordsieck yet and I'm not smart enough
   today to reverse-engineer his technique and try to make it work for the
   position-from-velocity case.  I'd rather go with a lower order
   adams-bashforth-moulton technique, since that's well documented and can do
   the position-from-velocity case. */
static void doalloc (Double **p, int size) {
  if (*p != 0) xfree (*p);
  *p = (Double *) xmalloc (size * sizeof (Double));
}

/* Copy the given state vector into the info. */
static void stateToInfo (const Double *const state,
                         BrennerMainInfo *const info, Double step) {
  const int num_atms = info->num_atms;
  int i;
  for (i = firstIndex (toState (info)); i < num_atms; i++) {
    const int pos = i * 6;
    BrenAtom *const a = &info->atm_num[i];
    a->coord.x = state [pos];
    a->coord.y = state [pos+1];
    a->coord.z = state [pos+2];
    a->velocity.x = state [pos+3] * step;
    a->velocity.y = state [pos+4] * step;
    a->velocity.z = state [pos+5] * step;
  }
}

static void eval (const Double *const state,
                  Double *const dstate,
                  BrennerMainInfo *const info, Double step) {
  BrenAtom *atm_num = info->atm_num;
  int i;
  const int num_atms = info->num_atms;
  stateToInfo (state, info, step);
#ifndef NDEBUG
  printf ("eval before find_atom_energies posx %f velx %f\n",
          info->atm_num[spy].coord.x, info->atm_num[spy].velocity.x);
#endif
  find_atom_energies(info);
#ifndef NDEBUG
  printf ("eval before find_atom_energies velx %f forcex %f\n",
          info->atm_num[spy].velocity.x, info->atm_num[spy].force.x);
#endif
  /* Copy the result into dstate. */
  {
    const int num_atms = info->num_atms;
    int i;
    for (i = firstIndex (toState (info)); i < num_atms; i++) {
      int pos = i * 6;
      BrenAtom *const a = &info->atm_num [i];
      const Double mass = info->atm_num[i].mass * ECONV;
      dstate [pos] = a->velocity.x/step;
      dstate [pos+1] = a->velocity.y/step;
      dstate [pos+2] = a->velocity.z/step;
      dstate [pos+3] = a->force.x / mass;
      dstate [pos+4] = a->force.y / mass;
      dstate [pos+5] = a->force.z / mass;
    }
  }
}

static void rungekutta (BrennerMainInfo *info) {
  Double m_step = info->timestep;
  const int m_startedSize = info->num_atms * 6;
  static int allocated_num_atms = -1;
  /* This closely follows fungimol/Numanal/RungeKutta.cpp.  It's my name on the
   copyright, so there's no issue about whether I can do this.  This copy of
   this source is contributed to Zyvex so it mixes smoothly with the other
   code in this file.  Tim Freeman 9/18/2000.  
   FIXME I haven't thought through yet when to use single precision here, so
   they will all be double for now.  Tim Freeman */
  static Double *state = 0;
  static Double *m_k1 = 0;
  static Double *m_k2 = 0;
  static Double *m_k3 = 0;
  static Double *m_k4 = 0;
  static Double *m_tmp = 0;
  int i;
  if (allocated_num_atms != m_startedSize) {
    allocated_num_atms = m_startedSize;
    doalloc (&state, m_startedSize);
    doalloc (&m_k1, m_startedSize);
    doalloc (&m_k2, m_startedSize);
    doalloc (&m_k3, m_startedSize);
    doalloc (&m_k4, m_startedSize);
    doalloc (&m_tmp, m_startedSize);
  }
  for (i = 0; i < m_startedSize; i++) {
    const int pos = i * 6;
    const BrenAtom *const a = &info->atm_num[i];
    state [pos] = a->coord.x;
    state [pos+1] = a->coord.y;
    state [pos+2] = a->coord.z;
    state [pos+3] = a->velocity.x/m_step;
    state [pos+4] = a->velocity.y/m_step;
    state [pos+5] = a->velocity.z/m_step;
  }
  eval (state, m_k1, info, m_step);
  for (i = 0; i < m_startedSize; i++) {
	m_tmp [i] = state [i] + m_step * m_k1 [i] / 2;
  }
  eval (m_tmp, m_k2, info, m_step);
  for (i = 0; i < m_startedSize; i++) {
	m_tmp [i] = state [i] + m_step * m_k2 [i] / 2;
  }
  eval (m_tmp, m_k3, info, m_step);
  for (i = 0; i < m_startedSize; i++) {
	m_tmp [i] = state [i] + m_step * m_k3 [i];
  }
  eval (m_tmp, m_k4, info, m_step);
  for (i = 0; i < m_startedSize; i++) {
	state [i] = state [i] +
	  m_step * (m_k1 [i] + 2 * m_k2 [i] + 2 * m_k3 [i] + m_k4 [i]) / 6;
  }
  stateToInfo (state, info, m_step);
}
#endif

void
bren_1_step(BrennerMainInfo *info, int kflag)
{
#ifdef RUNGEKUTTA
  /* We only do plain physics when we're doing Runge-Kutta.  The purpose of the
   Runge-Kutta code is to create numbers that can compare directly with what
   is generated by Fungimol.  I could imagine doing it the other way,
   changing Fungimol's numerical integration to match what is used here.  The
   problem with that plan is that I don't know how to do our third-order
   predictor-corrector in the case when the laws of physics are something
   other than F = MA, such as is the case for Fungimol's DesignAtom's.  Tim
   Freeman 18 Sep 2000. */
  rungekutta (info);
#else  
  third_order_predictor(info);
  find_atom_energies(info);
  
  if(kflag == -1)
    apply_gleq_thermostat(info);
  else if(kflag == 2) {
    int i;
    for(i = firstIndex (toState (info)); i < info->num_atms; ++i)
      info->atm_num[i].velocity.x =
        info->atm_num[i].velocity.y =
        info->atm_num[i].velocity.z = 0.0;
  } else if(kflag == 3)
    apply_hoov_thermostat(info);
  else if (kflag != 4) {
    if(info->num_thermostated_atoms > 0)
      apply_thermostat(info);
  }
  third_order_corrector(info);
  info->lchk = choose_lj(toState (info));
  if(info->zero_com_velocity)
    remove_com_velocity(info);
#ifndef INFINITE_CUBE
  if(info->volume_scale_dir != -1)
    vscale(info);
#endif
#endif
}

void
update_info(BrennerMainInfo *info)
{
}

BrennerMainInfo *
alloc_bren(int *kflag)
{
    BrennerMainInfo *info = (BrennerMainInfo*)malloc(sizeof(BrennerMainInfo));
    memset(info, 0, sizeof(*info));
    info->ljNeighborState = newNeighborState ();
    info->caNeighborState = newNeighborState ();
    info->temperature = 300;
    info->num_atms = 0;
    info->steps = 100000;
    info->timestep = 0.5;
    info->starttime = 0;
    info->lchk = 1;
    info->num_thermostated_atoms = 0;
    info->head[0] = 0;
    info->volume_scale_dir = -1;
    info->min_range_tree = 100;
    info->min_lj_range_tree = 100;
    initkt(info, kflag);
    init_c();
    PARAM(ktmax, info->RLL);
    return info;
}

void
init_bren(BrennerMainInfo *info, int kflag, int tight_binding)
{
	int i;
	int nonzero_vel = 0;
	for (i=firstIndex (toState (info)); i < info->num_atms; i++)
	{
		if(info->atm_num[i].thermostated)
		  ++info->num_thermostated_atoms;
		if(info->atm_num[i].velocity.x != 0 || 
		   info->atm_num[i].velocity.y != 0 ||
		   info->atm_num[i].velocity.z != 0)
		  nonzero_vel = 1;
	}
	if(!nonzero_vel) initialize_atoms(info, kflag);
	info->dnla=3*info->timestep*info->timestep*info->num_thermostated_atoms;
	if(kflag == 6)
	  minimize(info);
}
/* #endif for FUNGIMOL. */
#endif
