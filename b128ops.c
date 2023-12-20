/* Routines for arithmetic/logical operations on binary matrices 
   represented with a bitarray for each matrix row.
    
   These routines are a drop-in replacements (in k2comp & k2mult) 
   for those implementing k2-matrices. To ensure complete compatibility  
   their prototypes have some oddities:
     1) the variable asize is not used in the body of the function   
     2) the variable size usually contains the same vaue as in m->size
        except in matrix-creating functions mwrite_to_bbm() msave_to_file()

   Bit array representing rows are defined as arrays of uint128_t
   in future one could test the use of uint64_t and compare the speed 

   Matrix dimensions must be between 1 and 2^30

   Inside the b128mat_t the matrix size is represented as a size_t
   to ensure that when multiplied by a row/column index the result is 
   computed correctly using 64 bit arithmetic 


   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <string.h>
#include "b128.h"
#include "bbm.h"


static void quit(const char *msg, int line, char *file);
static void b128_free(b128mat_t *m);
static void b128_init(int size, b128mat_t *a);


// write the content of a b128 matrix :a to the bbm matrix :m 
// of size msize*msize. It is assumed m was already correctly initialized and allocated
void mwrite_to_bbm(uint8_t *m, int msize, int asize, const b128mat_t *a)
{
  (void) asize;
  assert(a!=NULL);
  if(msize!=a->size) quit("mwrite_to_bbm: matrix size mismatch", __LINE__,__FILE__);
  byte_to_bbm(m,msize,0,0,msize,2); // fill m with illegal value 2
  uint128_t one = 1;
  for(size_t i=0;i<msize;i++) {  // row index 
    size_t offset = i*a->colb;   // start of row i in a->b
    for(size_t j=0;j<msize;j++) 
      if(a->b[offset + j/128] & (one << (j%128)))
        m[i*msize+j] = 1;
      else 
        m[i*msize+j] = 0;
  }
}


// compress the matrix *m of size msize into the b128mat_t structure *a 
// m should be an  array of size msize*msize 
// the old content of :a is lost
// fro compatibilty with k2mat return the size of the b128 matrix (ie msize)
int mread_from_bbm(uint8_t *m, int msize, b128mat_t *a)
{
  assert(a!=NULL && m!=NULL);
  b128_free(a); // free previous content of a
  if(msize<=0) quit("mread_from_bbm: illegal matrix size",__LINE__,__FILE__);;
  if(msize>(1<<30)) quit("mread_from_bbm: matrix too large",__LINE__,__FILE__);
  b128_init(msize,a);
  bzero(a->b,a->size*a->colb*sizeof(uint128_t));
  uint128_t one = 1;
  for(size_t i=0;i<msize;i++) {
    size_t offset = i*a->colb;   // start of row i in a->b
    for(size_t j=0;j<msize;j++)
      if(m[i*msize+j]==1)
        a->b[offset + j/128] |= (one << (j%128));
      else
        assert(m[i*msize+j]==0);
  }
  return msize; // integer returned for compatibilty with k2mat
}

// write to :file statistics for b128_mat :a with an arbitrary :name as identifier
void mshow_stats(size_t size, int asize, const b128mat_t *a, const char *mname,FILE *file) {
  (void) asize;
  if(size!=a->size) quit("mshow_stats: size mismatch",__LINE__,__FILE__);
  fprintf(stderr,"%s -- matrix size: %zd, block size: %zd, # column blocks %d\n",
          mname,size,8*sizeof(*(a->b)), a->colb);  
}

// main entry point for matrix equality
// check if two b128 compressed matrices :a and :b are equal
// if a==b return -1
// if a!=b return the row index>=0 containing the first difference
int mequals(int size, const b128mat_t *a, const b128mat_t *b)
{
  (void) size;
  assert(a!=NULL && b!=NULL);
  assert(a->size!=b->size);
  for(size_t i=0; i<a->size*a->colb; i++)
    if(a->b[i]!=b->b[i]) return i/a->colb; // return row number of first difference
  return -1;
}



// matrix addition
// sum size x size b128 compressed matrices :a and :b storing
// the result in c, old content of c is discarded
// :a and :b must be of the same size
// Not tested!
void msum(int asize, const b128mat_t *a, const b128mat_t *b, b128mat_t *c)
{
  (void) asize;
  assert(a!=NULL && b!=NULL && c!=NULL);
  if(a->size != b->size)
    quit("msum: matrix size mismatch",__LINE__,__FILE__);

  // init and clean c
  b128_free(c);
  b128_init(a->size,c);
  assert(c->colb==a->colb && c->colb == b->colb);
  for(size_t i=0; i<a->size*a->colb; i++)
    c->b[i] = a->b[i] | b->b[i]; 
  return;
}

// main entry point for matrix multiplication. 
// multiply size x size b128 compressed matrices :a and :b storing
// the result in c, old content of c is discarded
// :a and :b must be of the same size
void mmult(int asize, const b128mat_t *a, const b128mat_t *b, b128mat_t *c)
{
  (void) asize;
  assert(a!=NULL && b!=NULL && c!=NULL);
  if(a->size != b->size)
    quit("mmult: matrix size mismatch",__LINE__,__FILE__);

  // init and clean c
  b128_free(c);
  b128_init(a->size,c);
  assert(c->colb==a->colb && c->colb == b->colb);
  bzero(c->b,c->size*c->colb*sizeof(uint128_t));
  uint128_t one = 1;
  for(size_t i=0; i<a->size; i++)  // row index for a 
    for(size_t k=0; k<a->size; k++) { // col index for a 
      size_t pos = i*a->colb+k/128;  
      if(a->b[pos] & (one << (k%128) )) { // a[i][k]==1
        // operates on whole rows:  c[i,:] |= b[k,:]
        for(int j=0; j<a->colb; j++)  // iterate over all row blocks 
          c->b[i*a->colb + j] |= b->b[k*a->colb +j];
      }
    }
}

// save the matrix :a to file :filename
// the format is
// 1) the size of the matrix (an int > 0)
// 2) the bit array (an uint128_t array of size a->size*a->colb)
void msave_to_file(int size, int asize, const b128mat_t *a, const char *filename)
{
  (void) asize;
  assert(a!=NULL);
  if(size!=a->size) quit("msave_to_file: matrix size mismatch", __LINE__,__FILE__);
  FILE *f = fopen(filename,"w");
  if(f==NULL) quit("msave_to_file: cannot open file", __LINE__,__FILE__);
  size_t e = fwrite(&size,sizeof(int),1,f);
  if(e!=1) quit("msave_to_file: cannot write size",__LINE__,__FILE__);
  size_t tot128 = a->size*a->colb;
  if(tot128 <=0) quit("msave_to_file: illegal bit array size",__LINE__,__FILE__);
  e = fwrite(a->b,sizeof(uint128_t),tot128,f);
  if(e!=tot128) quit("msave_to_file: cannot write bit array",__LINE__,__FILE__);
  fclose(f);
}


// load a b128 matrix stored in file :filename into the b128mat_t structure :a
// return the actual size of the matrix, see msave_to_file() for the file format
// :a old content is discarded
int mload_from_file(int *asize, b128mat_t *a, const char *filename)
{
  (void) *asize;
  assert(a!=NULL);
  b128_free(a);
  FILE *f = fopen(filename,"r");
  if(f==NULL) quit("mload_from_file: cannot open file", __LINE__,__FILE__);
  int size;
  size_t e = fread(&size,sizeof(int),1,f);
  if(e!=1) quit("mload_from_file: cannot read matrix size",__LINE__,__FILE__);
  if(size<1) quit("mload_from_file: matrix size smaller than 1, wrong format?",__LINE__,__FILE__);
  b128_init(size,a);
  size_t tot128 = a->size*a->colb;
  e = fread(a->b,sizeof(uint128_t),tot128,f);
  if(e!=tot128) quit("mload_from_file: cannot read bit array",__LINE__,__FILE__);
  fclose(f);
  return size;
}

// free memory used by m and init 
void matrix_free(b128mat_t *m) {
  b128_free(m);
}

// make a read-only copy of a matrix without allocating new memory
// old content of :c is discarded
void mmake_pointer(const b128mat_t *a, b128mat_t *c) {
  assert(a!=NULL && c!=NULL);
  b128_free(c);
  *c = *a; // copy all fields
  c->read_only = true;   // c is read only

}

// do nothing, added for compatibility with k2mat
void minimat_init(int m)
{
  (void) m;
}

// ----------- static auxiliary functions ------------

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

// free mem used by :m, old content is lost but *m still reusable if needed 
void b128_free(b128mat_t *m)
{
  if(m->read_only) // read only matrices are pointers to other matrices
    quit("Illegal operation: freeing a read only b128-matrix",__LINE__,__FILE__); 
  if(m->b!=NULL) free(m->b);
  m->b=NULL;
  m->size = m->colb = 0;
}

// make :a a b128 matrix of a given :size by initializing its fields
void b128_init(int size, b128mat_t *a)
{
  assert(a!=NULL);
  if(a->b!=NULL) quit("b128_init: initalizinga non-empty matrix",__LINE__,__FILE__);
  if(size<=0) quit("b128_init: illegal matrix size",__LINE__,__FILE__);
  a->size = size;
  a->colb = (size+127)/128;
  a->b = malloc(a->size*a->colb*sizeof(uint128_t));
  if(a->b==NULL) quit("b128_init: malloc failed",__LINE__,__FILE__);
  return;
}


