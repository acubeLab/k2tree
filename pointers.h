#ifndef _PAUXDEFS_H
#define _PAUXDEFS_H 1

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "vu64.h"

// type to represent the position of a node in a k2_tree
typedef uint32_t k2node_index_t;
// type to represent the position of a subtree info in k2->subtinfo 
typedef uint32_t k2sub_index_t; 

typedef struct pointers_t {
  // pointers info
  size_t p_size; // # of pointers
  k2node_index_t* nodep; // node pointers
  k2sub_index_t* subp;  // subtree pointers
} pointers_t;

void pointers_init(pointers_t **ps);
void pointers_copyinfo(pointers_t *ps, vu64_t* v);
void pointers_write_to_file(pointers_t *ps, const char* filename);
void pointers_load_from_file(pointers_t *ps, const char* filename);
void pointers_free(pointers_t* ps);
uint64_t pointers_size_in_bits(pointers_t* ps);

#endif
