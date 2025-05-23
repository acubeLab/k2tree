#ifndef _PAUXDEFS_H
#define _PAUXDEFS_H 1

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "vu64.h"

// type to represent the position of a node in a k2_tree
typedef uint64_t k2pointer_t;

typedef struct pointers_t {
  // pointers info
  k2pointer_t* nodep; // node pointers
  size_t     size; // # of pointers
} pointers_t;

pointers_t *pointers_init(vu64_t* v);
void pointers_write_to_file(pointers_t *ps, const char* filename);
pointers_t *pointers_load_from_file(const char* filename);
void pointers_free(pointers_t* ps);
uint64_t pointers_size_in_bits(pointers_t* ps);

#endif
