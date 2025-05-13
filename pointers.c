#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "pointers.h" 

static void quit(const char *msg, int line, char *file);

void pointers_init(pointers_t **ps) {
  *ps = malloc(sizeof(pointers_t));
  if(*ps == NULL) quit("malloc failed",__LINE__,__FILE__);
}

void pointers_copyinfo(pointers_t *ps, vu64_t* v) {
  assert(v->n > 0);
  ps->p_size = v->n;
  for(size_t i = 0; i < ps->p_size; i++) ps->p[i] = v->v[i];
}

void pointers_write_to_file(pointers_t *ps, const char* filename) {  
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening file", __LINE__, __FILE__);

  size_t check = fwrite(&(ps->p_size), sizeof(uint32_t), 1, file);
  if(check != sizeof(uint32_t)) quit("error writing the amount of pointers", __LINE__, __FILE__);
  check = fwrite(ps->p, sizeof(uint32_t), ps->p_size, file);
  if(check != ps->p_size) quit("error writing the pointers", __LINE__, __FILE__);
}

void poitners_load_from_file(pointers_t *ps, const char* filename) {  
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening file", __LINE__, __FILE__);

  size_t check = fread(&(ps->p_size), sizeof(uint32_t), 1, file);
  if(check != sizeof(uint32_t)) quit("error reading the amount of pointers", __LINE__, __FILE__);
  check = fwrite(ps->p, sizeof(uint32_t), ps->p_size, file);
  if(check != ps->p_size) quit("error reading the pointers", __LINE__, __FILE__);
}

void pointers_free(pointers_t* ps) {
  assert(ps != NULL);
  assert(ps->p != NULL);
  free(ps->p);
  free(ps);
}

uint64_t pointers_size_in_bits(pointers_t* ps) {
  if(ps == NULL) return 0;
  if(ps->p == NULL) return sizeof(ps->p_size);
  return sizeof(ps->p_size) + sizeof(ps->p[0]) * ps->p_size;
}

uint64_t pointers_size_in_bits(pointers_t* ps);
// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
