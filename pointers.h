#ifndef _PAUXDEFS_H
#define _PAUXDEFS_H 1

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "vu64.h"

typedef struct pointers_t {
  // pointers info
  uint32_t p_size; // amount of pointers
  uint32_t* p; // pointers
} pointers_t;

void pointers_init(pointers_t **ps);
void pointers_copyinfo(pointers_t *ps, vu64_t* v);
void pointers_write_to_file(pointers_t *ps, const char* filename);
void poitners_load_from_file(pointers_t *ps, const char* filename);
void pointers_free(pointers_t* ps);
uint64_t pointers_size_in_bits(pointers_t* ps);

#endif
