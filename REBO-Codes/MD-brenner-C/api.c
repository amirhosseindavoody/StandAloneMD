/*
 This code was released into the public domain by Tim Freeman on 11 Aug 2000.
*/

#ifndef __ANYAPI_H__
#include "anyapi.h"
#endif

void setAllPrevToCurrent (struct State *s) {
  int last = pastLastIndex (s);
  int i;
  for (i = firstIndex (s); i < last; i++) {
    setPrevToCurrent (i, s);
  }
}
