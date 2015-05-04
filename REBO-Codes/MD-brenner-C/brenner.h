#ifndef __brenner_h__
#define __brenner_h__

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

#ifndef __BRENNER_H__
#define __BRENNER_H__

#ifndef __myassert_h__
#include "myassert.h"
#endif

#include "vector.h"
#include <stdio.h>

/* NTAB is the number of points to use in interpolation tables.
   FIXME If we made this smaller we'd get more cache hits and maybe not lose
   significant accuracy.  Maybe worth trying someday. */
#define NTAB 10000

#ifndef FUNGIMOL
#define MAX_ATOMS 15000

#define MAXNBORS 100
#endif

#define HYDROGEN 1
/*#define CARBON 12 */
#define CARBON 6
#define SILICON 14
#define GERMANIUM 32

#ifndef __KT_H__
#include "kt.h"
#endif

#ifndef __BOOL_H__
#include "bool.h"
#endif

/* Stop if something is a NAN. */
static inline Double doubleCheckNAN (Double d) {
#ifndef NDEBUG
  assert ((d < 0.0) || (d >= 0.0));
#endif
  return d;
}

static inline Float floatCheckNAN (Float d) {
#ifndef NDEBUG
  assert ((d < 0.0) || (d >= 0.0));
#endif
  return d;
}

#ifndef FUNGIMOL
/* An atom has coordinates, a velocity vector, an acceleration vector,
   a number of bonds to H and C, a possibility of being in one of several
   types of conjugated systems, an energy dependent on its bonded neighbors'
   location, and a list of neighbors. If it is not movable, don't move it. */

struct atm {
  dvector coord; /* Position, in angstroms. */
  /* The velocity field is really the velocity in angstroms per femtosecond
     times the timestep, which is measured in femtoseconds. */
  vector velocity;
  /* accel is perhaps the acceleration in angstroms per square femtosecond,
     times the timestep squared, divided by 2. */
  vector accel;
  /* dx3dt3 is perhaps the jerk in angstroms per cubic femtosecond,
     times the timestep cubed, divided by 3. */
  vector dx3dt3;
  /* If you take the force that's computed, divide by the mass field, and
     multiply by ECONV, then you get acceleration measured in angstroms per
     square femtosecond. */
  vector force;
  int number;			/* id used as index into array of atoms */
  int movable;
  /* Type is the atomic number. */
  int type;
  BOOL thermostated;
  /* ktype is the index into the kt2 array, for computing the lennard-jones
     potential.  The first four elements are wired into the code as:
     Carbon 1
     Hydrogen 2
     Silicon 3
     Germanium 4
     in file bren2.c function initkt and then the rest are allocated
     differently from run to run based on what is in input.d.
  */
  int ktype;
  /* The unit for the mass field is the weight of a hydrogen atom. */
  Float mass; 
  vector prev_coord;
};
typedef struct atm BrenAtom;

#endif

#ifndef __NEIGHBORSTATE_H__
#include "NeighborState.h"
#endif

#ifndef __INTPAIR_H__
#include "IntPair.h"
#endif

/* Get coefficients for H and F and torsion. */
void PARAM(int ktmax, Float RLL);

#ifndef FUNGIMOL
typedef struct
{
  BrenAtom atm_num[MAX_ATOMS];
  int num_atms;
  int nxmol;
#ifndef INFINITE_CUBE
  Float cube[3];
#endif
  int noa[NTYPES];	/* number of atoms by type */
  /* This field is set to 3 * timestep * timestep * number of thermostatted
     atoms in init_bren. */ 
  Float dnla;
  Float temperature; /* kelvin */
  Float PSEED; /* Random # seed */
  Float RLL; /* added radius for neighbor threshold */
  Float timestep; /* time in fs between steps */ 
  Float starttime;
  int steps;
  int num_thermostated_atoms;
  /* lchk is computed by choose_lj in ljguts.c.  It is set to 2 if nobody moved
     more than 0.5 * RLL from their previous coordinates.  It is set to 1 if we
     need to recompute the neighbor list, and in that case the previous
     coordinates are set to the current coordinates.

     Unless we're adding 0.5 * RLL of slop to our neighboring, this implies
     that our neighboring is an approximation.  Tim Freeman  5 Sep 2000.
  */
  int lchk;
  Float system_energy;
  int zero_com_velocity; /* if no rigid atoms, reset the velocity of molecule's
			    center of mass to zero after each step */
  int volume_scale_dir; /* -1 for no scaling, or 0, 1, or 2 to scale
			 the volume along the x, y, or z axis */
  int squeeze;           /* set to true to speed up neighbor code when
			    there is a large number of non-movable atoms */
  int min_lj_range_tree; /* min #atoms for using range tree for lj neighbors */
  int min_range_tree; /* min #atoms for using range tree for non-lj neighbors */
  /* head is read from the first line of coord.d, and it's written whenever we
     write coord.d and as the first line of each block of xmol.d */
  char head[512]; 
  const char *param_file_dir; /* Where *.d files are.  Must end in /
                                 if not empty. Works only with the
                                 Fortran library at the moment. */
  /* We carry ljNeighborState here instead of making it local to ljguts.c
     because eventually I want the force field computation
     in ljguts.c to be unconcerned with how the neighbors are computed.  
     Fungimol's neighboring code is incremental from one timestep to the next,
     so we'll want to carry around information about neighboring in the state.
     It's best to only carry around one thing, so I'm sticking everything that
     might have to persist into the BrennerMainInfo structure.  Tim Freeman 17
     Aug 2000. */
  /* ljNeighborState and caNeighborState are initialized to new NeighborState's
   by alloc_bren. */
  struct NeighborState *ljNeighborState;
  struct NeighborState *caNeighborState;
} BrennerMainInfo;

void minimize(BrennerMainInfo *info);
BrennerMainInfo *alloc_bren(int *kflag);
void init_bren(BrennerMainInfo *info, int kflag, int tight_binding);
void bren_1_step(BrennerMainInfo *info, int kflag);
void update_info(BrennerMainInfo *info);
int set_cell_threshold(int n);
#endif

void my_exit(int);

#ifndef NDEBUG
/* The atom number we want to spy on.  This is the position in atm_num, which
   is one less than the atom number appearing in xmol.d.  Set this to -1 if you
   want less chatter in debug builds. */
#define spy -1
#endif

#ifndef BRENFORT_H
#ifndef __BRENC_H__
#include "brenc.h"
#endif
int out_of_box(double);
#endif

#endif

#endif
