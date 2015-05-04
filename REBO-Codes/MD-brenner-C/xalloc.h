/*
 This code was released into the public domain by Tim Freeman on 20 Aug 2000.

 Memory allocation with error checking.
*/

#ifndef __XALLOC_H__
#define __XALLOC_H__
/* Next one is for size_t and free. */
#include <stdlib.h>

#ifdef FUNGIMOL
#define FUNGIMOL_MEMORY 1
#endif

#ifdef FUNGIMOL_MEMORY

#ifndef __MemoryUtil_h__
#include "MemoryUtil.h"
#endif

#define xrealloc(ptr,size) realloc (ptr, size)
#define xmalloc(size) malloc(size)
/* The next one is for consistency only.  By calling xfree we make it
   impossible for people to compile without including xalloc.h, and then it's
   impossible to have subtle failures when FUNGIMOL_MEMORY is turned on.  The
   subtle failure is that if we use Fungimol's malloc and the system's free, we
   seg fault when we free. */
#define xfree(ptr) free(ptr)
#else
/* Stop if reallocation failed. */
void *xrealloc (void *ptr, size_t size);
/* Stop if allocation failed. */
void *xmalloc (size_t size);
/* The next one is for consistency only.  By calling xfree we make it
   impossible for people to compile without including xalloc.h, and then it's
   impossible to have subtle failures when FUNGIMOL_MEMORY is turned on. */
#define xfree(ptr) free(ptr)
#endif

#endif
