#ifndef _PAUXDEFS_H
#define _PAUXDEFS_H 1

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "vu64.h"

// #define SIMPLEBACKPOINTERS  // use simple backpointers, not the full ones

// type to represent the position of a node in a k2_tree
#ifdef SIMPLEBACKPOINTERS
typedef uint32_t k2pointer_t;
#define MAXPOINTER UINT32_MAX  // maximum value for a pointer: largets 40-bit value
#else
typedef uint64_t k2pointer_t;
#define MAXPOINTER TSIZEMASK   // maximum value for a pointer: largets 40-bit value
#endif

typedef struct pointers_t {
  // pointers info
  k2pointer_t *nodep;   // node pointers
  size_t        size;   // # of pointers
  uint32_t   *sorted;   // sorted order of pointers by destination, used only for construction
  size_t        sidx;   // index in the sorted array
} pointers_t;

pointers_t *pointers_init(vu64_t* v);
void pointers_write_to_file(pointers_t *ps, const char* filename);
pointers_t *pointers_load_from_file(const char* filename);
void pointers_free(pointers_t* ps);
uint64_t pointers_size_in_bits(pointers_t* ps);
void pointers_sort(pointers_t* ps);

#endif
