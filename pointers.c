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
  assert(ps != NULL);
  assert(v != NULL);
  assert(v->v != NULL);
  assert(v->n > 0);
  ps->p_size =  v->n;
  ps->nodep = malloc(sizeof(k2node_index_t) * ps->p_size);
  if(ps->nodep == NULL) quit("malloc failed",__LINE__,__FILE__);
  for(size_t i = 0; i < ps->p_size; i++) ps->nodep[i] = (k2node_index_t) v->v[i];
}

// write the node pointers to file
void pointers_write_to_file(pointers_t *ps, const char* filename) {  
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening file", __LINE__, __FILE__);

  size_t check = fwrite(ps->nodep, sizeof(k2node_index_t), ps->p_size, file);
  if(check != ps->p_size) quit("error writing the pointers", __LINE__, __FILE__);
  fclose(file);
}


//\\!!!! this function should allocate the pointer array
void pointers_load_from_file(pointers_t *ps, const char* filename) {  
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening pointers file", __LINE__, __FILE__);

  size_t check = fread(&(ps->p_size), sizeof(uint32_t), 1, file);
  if(check != 1) quit("error reading the number of pointers", __LINE__, __FILE__);
  check = fread(ps->nodep, sizeof(uint32_t), ps->p_size, file);
  if(check != ps->p_size) quit("error reading the pointers", __LINE__, __FILE__);
  fclose(file);
}

void pointers_free(pointers_t* ps) {
  assert(ps != NULL);
  assert(ps->nodep != NULL);
  free(ps->nodep);
  free(ps);
}

uint64_t pointers_size_in_bits(pointers_t* ps) {
  if(ps == NULL) return 0;
  // there are no pointers
  if(ps->nodep == NULL) return sizeof(*ps);
  return (sizeof(*ps) + sizeof(ps->nodep[0]) * ps->p_size) * 8;
}

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
