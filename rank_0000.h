#ifndef _RANK0000DEFS_H
#define _RANK0000DEFS_H 1

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct rank_0000_t {
  // rank support 
  uint32_t r_size; // size of rank array
  uint32_t block_size; // block size
  uint32_t *r; // rank array, that is actually just a prefix sum
} rank_0000_t;

// the void is a k2mat_t, is just to avoid compilation conflicts
void rank_init(rank_0000_t **r, uint32_t block_size, void *a);
uint32_t rank_rank(rank_0000_t* r, void* a, uint32_t i);

void rank_write_to_file(rank_0000_t *r, const char* filename);
void rank_load_from_file(rank_0000_t *r, const char* filename);

void rank_free(rank_0000_t* r);
uint64_t rank_size_in_bits(rank_0000_t* r);

#endif
