#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "k2.h" 

static void quit(const char *msg, int line, char *file);

// init a pointers structure using data from the v array
pointers_t *pointers_init(vu64_t* v) {
  assert(v != NULL);
  assert(v->v != NULL);
  assert(v->n > 0);
  pointers_t *ps = malloc(sizeof(pointers_t));
  if(ps == NULL) quit("malloc failed",__LINE__,__FILE__);
  ps->size =  v->n;
  ps->nodep = malloc(sizeof(k2pointer_t) * ps->size);
  if(ps->nodep == NULL) quit("malloc failed",__LINE__,__FILE__);
  for(size_t i = 0; i < ps->size; i++) {
    if(v->v[i] > MAXPOINTER) {
      fprintf(stderr, "error: pointer value %lu exceeds maximum value %lu\n",
              (unsigned long)v->v[i], (unsigned long)MAXPOINTER);
      exit(EXIT_FAILURE);
    }
    ps->nodep[i] = (k2pointer_t) v->v[i];
  }
  ps->sorted = NULL; // no sorted order yet
  ps->sidx = 0; // no sorted order yet
  return ps;
}


// write the node pointers info to file
void pointers_write_to_file(pointers_t *ps, const char* filename) { 
  assert(ps!=NULL); 
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening file", __LINE__, __FILE__);

  size_t check = fwrite(ps->nodep, sizeof(k2pointer_t), ps->size, file);
  if(check != ps->size) quit("error writing pointers to file", __LINE__, __FILE__);
  fclose(file);
}


pointers_t *pointers_load_from_file(const char* filename) {  
  FILE* file = fopen(filename, "rb");
  if(file == NULL) quit("error opening pointers file", __LINE__, __FILE__);
  // read the size of the file
  fseek(file, 0, SEEK_END);
  long size = ftell(file);
  if(size < 0) quit("error getting pointers file size", __LINE__, __FILE__);
  // check if the file is empty
  if(size == 0) quit("error: pointers file is empty", __LINE__, __FILE__);

  // check if the file size is a multiple of the size of a pointer
  if(size % sizeof(k2pointer_t) != 0) 
    quit("error: file size is not a multiple of pointer size", __LINE__, __FILE__);
  // allocate and read pointers
  pointers_t *ps = malloc(sizeof(pointers_t));
  if(ps == NULL) quit("malloc failed",__LINE__,__FILE__);
  ps->size = size / sizeof(k2pointer_t);
  // allocate memory for the pointers
  ps->nodep = malloc(sizeof(k2pointer_t) * ps->size);
  if(ps->nodep == NULL) quit("malloc failed", __LINE__, __FILE__);
  // read the pointers from the file
  rewind(file);
  size_t check = fread(ps->nodep, sizeof(k2pointer_t), ps->size, file);
  if(check != ps->size) quit("error reading the pointers", __LINE__, __FILE__);
  fclose(file);
  ps->sorted= NULL; // no sorted order used yet
  ps->sidx = 0; 
  return ps;
}

// free structure representing pointer 
void pointers_free(pointers_t* ps) {
  assert(ps != NULL);
  if(ps->nodep != NULL)  free(ps->nodep);
  if(ps->sorted != NULL) free(ps->sorted);
  ps->nodep = NULL;
  ps->size = 0;
  ps->sorted = NULL;
  free(ps);
}

// return space usage of pointers structure in bits
// do not include the size of sorted which is used only for construction
uint64_t pointers_size_in_bits(pointers_t* ps) {
  if(ps == NULL) return 0;
  // there are no pointers
  if(ps->nodep == NULL) return sizeof(*ps)*8;
  return (sizeof(*ps) + sizeof(k2pointer_t) * ps->size) * 8;
}

// compare function for qsort_r
static int pointers_cmp(const void *a, const void *b, void *arg) {
  k2pointer_t *p = (k2pointer_t *) arg;
  uint32_t i = *(uint32_t *) a;
  uint32_t j = *(uint32_t *) b;
  if(p[i] < p[j]) return -1;
  else if(p[i] > p[j]) return 1;
  return 0;
}

// initialze the field ps->sorted containing the order of the pointers
// when sorted according to their destination
// before the sorting clear the 24 higher bits of the pointers
void pointers_sort(pointers_t* ps) {
  assert(ps != NULL);
  assert(ps->nodep != NULL && ps->size > 0);
  if(ps->size >= UINT32_MAX) 
    quit("error: too many pointers", __LINE__, __FILE__);
  // clear the 24 higher bits of the pointers
  for(size_t i = 0; i < ps->size; i++) {
    ps->nodep[i] &= ((k2pointer_t) 1 << BITSxTSIZE) - 1; // clear higher bits
  }
  // create permutation array ordering pointers by destination
  ps->sorted = malloc(sizeof(uint32_t) * ps->size);
  if(ps->sorted == NULL) quit("malloc failed", __LINE__, __FILE__);
  for(uint32_t i = 0; i < ps->size; i++) {
    ps->sorted[i] = i;
  }
  qsort_r(ps->sorted, ps->size, sizeof(uint32_t), pointers_cmp, ps->nodep);
  #ifndef NDEBUG
  // check that the pointers are sorted
  for(uint32_t i = 1; i < ps->size; i++) {
    if(ps->nodep[ps->sorted[i]] < ps->nodep[ps->sorted[i-1]]) {
      fprintf(stderr, "error: pointers are not sorted\n");
      exit(EXIT_FAILURE);
    }
  }
  #endif
  ps->sidx = 0; // reset index in the sorted array
}


// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
