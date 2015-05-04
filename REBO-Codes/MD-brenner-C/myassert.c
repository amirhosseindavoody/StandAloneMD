/*
 This copy of this code was released into the public domain by Tim Freeman on
 4 Oct 2000.  Never mind that Fungimol has a very similar file.  It's my
 copyright so I can do that.
*/
#include "myassert.h"

/* For abort. */
#ifndef __stdlib_h__
#include <stdlib.h>
#define __stdlib_h__
#endif

/* For printf. */
#ifndef __stdio_h__
#include <stdio.h>
#define __stdio_h__
#endif

int myassertimplementation
(const char *string, const char *file, const int line,
 const char *pretty_function) {
  fprintf (stderr, "Assertion failed at line %d of %s in %s:\n%s\n",
           line, file, pretty_function, string);
  abort ();
  return 1;
}
