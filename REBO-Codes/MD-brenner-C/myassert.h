/*
 This copy of this code was released into the public domain by Tim Freeman on
 4 Oct 2000.  Never mind that Fungimol has a very similar file.  It's my
 copyright so I can do that.
*/

#ifndef __myassert_h__
#define __myassert_h__

#if 1
   int myassertimplementation
   (const char *string, const char *file, const int line,
    const char *pretty_function);
   /* GCC's assert in assert.h garbages the caller's stack frame when it
   triggers, making debugging difficult.
   If they included <assert.h>, get rid of the old assert macro. */
   #ifdef assert
   /* I think I spelled it "systems" instead of "system's" in the next line of
   code because my editor or the compiler or something got confused by the
   stray single quote. */
      #error You included the systems assert.h!
   #endif
   #undef assert
   // Next line should cause the external include guards not to
   // redundantly include <assert.h>  
   #define __assert_h__
   #ifdef NDEBUG
      #define assert(x) ((void) 0)
   #else
      #define assert(x) ((void) ((x)?0: \
		      myassertimplementation(#x, __FILE__, __LINE__, \
			       __PRETTY_FUNCTION__)))
   #endif
#else
   #ifndef __assert_h__
   #include <assert.h>
   #define __assert_h__
   #endif
#endif
#endif
