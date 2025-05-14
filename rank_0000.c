#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "k2.h"
#include "rank_0000.h" 

static void quit(const char *msg, int line, char *file);
node_t k2read_node__(const k2mat_t *m, size_t p); // maybe later I can find a better way

void rank_init(rank_0000_t **r, uint32_t block_size, void *a) {
  k2mat_t* a_ = (k2mat_t*) a;

  *r = malloc(sizeof(rank_0000_t));
  (*r)->r_size = (uint32_t) (a_->pos + block_size - 1) / block_size + 1;
  (*r)->block_size = block_size;
  (*r)->r = (uint32_t*) malloc(sizeof(uint32_t) * ((*r)->r_size));

  uint32_t sum = 0;
  for(size_t i = 0; i < a_->pos; i++) {
    if(i % block_size == 0) {// finish block
      (*r)->r[(i + block_size - 1) / block_size] = sum;
    }
    sum += k2read_node__(a_, i) == POINTER;
  }

  (*r)->r[(a_->pos + block_size - 1) / block_size] = sum;
}

uint32_t rank_rank(rank_0000_t* r, const void* a, uint32_t i) {
  k2mat_t* a_ = (k2mat_t*) a;
  assert(i <= a_->pos);
  uint32_t block = (uint32_t) i / r->block_size;

  assert(block < r->r_size);
  uint32_t ret = r->r[block];

  for(uint32_t to_read = block * r->block_size; to_read < i; to_read++) {
    ret += k2read_node__(a_, to_read) == 0;
  }
  return ret;
}

void rank_write_to_file(rank_0000_t *r, const char* filename) {
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening file", __LINE__, __FILE__);

  size_t check = fwrite(&(r->r_size), sizeof(uint32_t), 1, file);
  if(check != 1) quit("error writing the rank size", __LINE__, __FILE__);
  check = fwrite(&(r->block_size), sizeof(uint32_t), 1, file);
  if(check != 1) quit("error writing the block size", __LINE__, __FILE__);
  check = fwrite(r->r, sizeof(uint32_t), r->r_size, file);
  if(check != r->r_size) quit("error writing the prefix sum", __LINE__, __FILE__);
  fclose(file);
}

void rank_load_from_file(rank_0000_t *r, const char* filename) {
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening file", __LINE__, __FILE__);

  size_t check = fread(&(r->r_size), sizeof(uint32_t), 1, file);
  if(check != 1) quit("error reading the rank size", __LINE__, __FILE__);
  check = fread(&(r->block_size), sizeof(uint32_t), 1, file);
  if(check != 1) quit("error reading the block size", __LINE__, __FILE__);
  check = fwrite(r->r, sizeof(uint32_t), r->r_size, file);
  if(check != r->r_size) quit("error reading the pointers", __LINE__, __FILE__);
  fclose(file);
}

void rank_free(rank_0000_t* r) {
  assert(r != NULL);
  assert(r->r != NULL);
  free(r->r);
  free(r);
}

uint64_t rank_size_in_bits(rank_0000_t* r) {
  if(r == NULL) return 0;
  return (sizeof(r->r_size) + 
          sizeof(r->block_size) + 
          sizeof(r->r) + sizeof(r->r[0]) * r->r_size) * 8;
}

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

node_t k2read_node__(const k2mat_t *m, size_t p)
{
  p+=m->offset;
  assert(p<m->pos && p<m->lenb);
  if(p%2==0)
    return m->b[p/2] & 0xF;
  else 
    return  (m->b[p/2] >> 4) & 0xF;
}
