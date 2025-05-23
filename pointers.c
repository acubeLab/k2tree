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

// init a pointers structure usng data from the v array
pointers_t *pointers_init(vu64_t* v) {
  assert(v != NULL);
  assert(v->v != NULL);
  assert(v->n > 0);
  pointers_t *ps = malloc(sizeof(pointers_t));
  if(ps == NULL) quit("malloc failed",__LINE__,__FILE__);
  ps->size =  v->n;
  ps->nodep = malloc(sizeof(k2pointer_t) * ps->size);
  if(ps->nodep == NULL) quit("malloc failed",__LINE__,__FILE__);
  for(size_t i = 0; i < ps->size; i++) ps->nodep[i] = (k2pointer_t) v->v[i];
  return ps;
}


// write the node pointers info to file
void pointers_write_to_file(pointers_t *ps, const char* filename) {  
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening file", __LINE__, __FILE__);

  size_t check = fwrite(ps->nodep, sizeof(k2pointer_t), ps->size, file);
  if(check != ps->size) quit("error writing pointers to file", __LINE__, __FILE__);
  fclose(file);
}


pointers_t *pointers_load_from_file(const char* filename) {  
  FILE* file = fopen(filename, "w");
  if(file == NULL) quit("error opening pointers file", __LINE__, __FILE__);
  // read the size of the file
  fseek(file, 0, SEEK_END);
  long size = ftell(file);
  if(size < 0) quit("error getting file size", __LINE__, __FILE__);
  // check if the file is empty
  if(size == 0) {
    fclose(file);
    return NULL;
  }
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
  return ps;
}

// free structure reprsenting pointer 
void pointers_free(pointers_t* ps) {
  assert(ps != NULL);
  if(ps->nodep != NULL) {
    free(ps->nodep);
    ps->nodep = NULL;
    ps->size = 0;
  }
  free(ps);
}

uint64_t pointers_size_in_bits(pointers_t* ps) {
  if(ps == NULL) return 0;
  // there are no pointers
  if(ps->nodep == NULL) return sizeof(*ps)*8;
  return (sizeof(*ps) + sizeof(k2pointer_t) * ps->size) * 8;
}

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
