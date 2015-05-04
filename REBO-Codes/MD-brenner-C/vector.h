/* vector.c - Copyright (c) 1998 Zyvex LLC.
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
 *  This product includes software developed by Zyvex LLC.
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

#ifndef __VECTOR_H__
#define __VECTOR_H__

/* The next one is required for pow. */
#ifndef MATH_H
#include <math.h>
#define MATH_H
#endif

typedef struct { Float x,y,z; } vector;
typedef struct { Double x,y,z; } dvector;

/* Vectors used to require a subroutine call to use.  Now I inline them in this
   header file.  The duration of the timeruns changed from 4.90 seconds to
   4.88.  I had hoped for more dramatic results, but given that there were any
   results at all I'm unwilling to un-inline.  Inlining vectors made a much
   larger difference in Fungimol, I think because I used the analogous vector
   class more there.  Tim Freeman 29 Jul 2000. */

/* Add two vectors */
static inline vector plus( vector i, vector j )
{
    vector sum;
    sum.x = i.x + j.x;
    sum.y = i.y + j.y;
    sum.z = i.z + j.z;
    return sum;
}

static inline dvector dplus( dvector i, vector j )
{
    dvector sum;
    sum.x = i.x + j.x;
    sum.y = i.y + j.y;
    sum.z = i.z + j.z;
    return sum;
}

static inline vector minus ( vector i, vector j )
{
    vector sum;
    sum.x = i.x - j.x;
    sum.y = i.y - j.y;
    sum.z = i.z - j.z;
    return sum;
}

/* Note that subtracting two double-precision vectors gives you a
   single-precision vector.  This is because we use double-precision vectors
   for position, and I don't think we need the precision for the difference
   between two vectors. */
static inline vector ddminus ( dvector i, dvector j )
{
    vector sum;
    sum.x = i.x - j.x;
    sum.y = i.y - j.y;
    sum.z = i.z - j.z;
    return sum;
}

/* If we switch to C++, we don't have to keep straight the difference in
   meaning between plus and sum, and we wouldn't have to have the different
   spellings of minus and plus above.  But IMO it isn't worth switching to C++
   only for that.  Tim Freeman 30 Jul 2000. */
static inline vector sum( Float i, vector j )
{
    vector result;
    result.x = i + j.x;
    result.y = i + j.y;
    result.z = i + j.z;
    return result;
}

/* Find the product of a scalar and a vector. */
static inline vector product( Float i, vector j )
{
    vector vec;
    vec.x = i*j.x; 
    vec.y = i*j.y; 
    vec.z = i*j.z;
    return vec;
}

/* Raise a vector to a power elementwise. */

static inline vector power( Float i, vector j )
{
    vector vec;
    vec.x = pow(j.x, i);
    vec.y = pow(j.y, i);
    vec.z = pow(j.z, i);
    return vec;
}

/* Find the dot product of two vectors. */
static inline Float dot( vector i, vector j )
{
    return ( i.x*j.x + i.y*j.y + i.z*j.z );
}

/* Find the cross product of two vectors */
static inline vector cross( vector b, vector c)
{
    vector vec;
    vec.x = b.y*c.z - c.y*b.z;
    vec.y = b.x*c.z - c.x*b.z;
    vec.z = b.x*c.y - c.x*b.y;
    return vec;
}

#endif
