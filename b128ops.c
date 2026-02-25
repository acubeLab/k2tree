/* Routines for arithmetic/logical operations on binary matrices 
   represented with a bitarray for each matrix row.
    
   These routines are a drop-in replacements (in k2comp & k2mult) 
   for those implementing k2-matrices. To ensure complete compatibility  
   their prototypes have some oddities:
     1) the variable asize is not used in the body of the function   
     2) the variable size usually contains the same value as in m->size
        except in matrix-creating functions mwrite_to_bbm() msave_to_file()

   Bit array representing rows are defined as arrays of uint128_t
   in future one could test the use of uint64_t and compare the speed
   The bits outsize the matrix in the last 128-bit blocks are
   guaranteed to be zero.  

   Matrix dimensions must be between 1 and 2^30 (in practice much less
   since a size 2^30 matrix would take 2^57 bytes)  

   However, in the b128mat_t def the matrix size is represented as a size_t
   for compatibility with sparse representations and to ensure that when 
   the size is multiplied by a row/column index the result 
   is computed correctly using 64 bit arithmetic 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <errno.h>
#include <assert.h>
#include <string.h>
#include <inttypes.h>
#include "b128.h"
#include "bbm.h"

static void quit(const char *msg, int line, char *file);
static void b128_free(b128mat_t *m);
static void b128_init(size_t size, b128mat_t *a);


// read a matrix from the text file :iname (one entry per line)
// of size xsize and store it into the b128mat_t structure *a 
// the old content of :a is lost
// for compatibilty with k2mats return the size of the b128 matrix (ie msize)
// and store the same value to *msize
size_t mread_from_textfile(b128mat_t *a, char *iname, size_t xsize) {
  assert(a!=NULL && iname!=NULL);
  b128_free(a); // free previous content of a
  if(xsize<=0) quit("mread_from_textfile: illegal matrix size",__LINE__,__FILE__);;
  if(xsize>MaxMatrixSize) quit("mread_from_textfile: matrix too large",__LINE__,__FILE__);
  b128_init(xsize,a);
  bzero(a->b,a->size*a->colb*sizeof(uint128_t));
  uint128_t one = 1;
  // read nonzero entries from file
  FILE *f = fopen(iname,"rt");
  if(f==NULL) quit("mread_from_file: cannot open input file",__LINE__,__FILE__);
  int64_t i,j; size_t line=0;  
  while(true) {
    line++;
    int e = fscanf(f,"%" SCNd64 " %" SCNd64,&i,&j);
    if(e==EOF) break;
    // check input
    if(e!=2) {
      fprintf(stderr,"Invalid file content at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(i<0 || j<0) {
      fprintf(stderr,"Negative index at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(i>=xsize || j>=xsize) {
      fprintf(stderr,"Index larger than the assigned size at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    // write 1 in approprite bit position 
    a->b[i*a->colb + j/128] |= (one << (j%128));
  }
  fclose(f);
  a->size = xsize;
  return xsize;
}

// write the content of a b128 matrix :a to 
// text file in one entry per line format
void mwrite_to_textfile(const b128mat_t *a, char *outname)
{
  assert(outname!=NULL && a!=NULL);
  FILE *f = fopen(outname,"wt");
  if(f==NULL) quit("mwrite_to_file: cannot open output file",__LINE__,__FILE__);
  
  uint128_t one = 1;
  for(size_t i=0;i<a->size;i++) {  // row index 
    size_t offset = i*a->colb;   // start of row i in a->b
    for(size_t j=0;j<a->size;j++) 
      if(a->b[offset + j/128] & (one << (j%128)))
        fprintf(f,"%zu %zu\n",i,j);
  }
  fclose(f);
}



// write the content of a b128 matrix :a to the bbm matrix :m 
// of size msize*msize. It is assumed m was already correctly initialized and allocated
void mwrite_to_bbm(uint8_t *m, const b128mat_t *a)
{

  assert(a!=NULL);
  size_t msize = a->size;
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


// compress the bbm matrix *m of size msize into the b128mat_t structure *a 
// m should be an array of size msize*msize 
// the old content of :a is lost
// for compatibilty with k2mats return the size of the b128 matrix (ie msize)
size_t mread_from_bbm(uint8_t *m, size_t msize, b128mat_t *a)
{
  assert(a!=NULL && m!=NULL);
  b128_free(a); // free previous content of a
  if(msize<=0) quit("mread_from_bbm: illegal matrix size",__LINE__,__FILE__);;
  if(msize>MaxMatrixSize) quit("mread_from_bbm: matrix too large",__LINE__,__FILE__);
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
size_t mshow_stats(const b128mat_t *a, const char *mname,FILE *file) {
  fprintf(stderr,"%s:\n matrix size: %zu, block size: %zu, column blocks: %u\n",
          mname,a->size,8*sizeof(*(a->b)), a->colb);  
  return 0; // integer returned for compatibilty with k2mat (should be # nonzeros)
}

// main entry point for matrix equality
// check if two b128 compressed matrices :a and :b are equal
// if a==b return -1
// if a!=b return the row index>=0 containing the first difference
int mequals_plain(size_t size, const b128mat_t *a, const b128mat_t *b)
{
  (void) size;
  assert(a!=NULL && b!=NULL);
  assert(a->size==b->size);
  for(size_t i=0; i<a->size*a->colb; i++)
    if(a->b[i]!=b->b[i]) return (int) (i/a->colb); // return row number of first difference
  return -1;
}

// as above but return true if a and b are equal, false otherwise
bool mequals(const b128mat_t *a, const b128mat_t *b) {
  assert(a!=NULL && b!=NULL);
  if(a->size != b->size) return false; // cannot say
  return mequals_plain(0, a, b) < 0;
}

// add indentity matrix to a
void madd_identity(b128mat_t *a)
{
  uint128_t one = 1;
  for(size_t i=0; i<a->size; i++) 
    a->b[i*a->colb + i/128] |= (one << (i%128));
}


// creates a size x size zero matrix
b128mat_t mat_zero(b128mat_t *b) {
  b128mat_t a = B128MAT_INITIALIZER;
  b128_init(b->size,&a);
  bzero(a.b,a.size*a.colb*sizeof(uint128_t));
  return a;
}


// creates a size x size identity matrix
b128mat_t mat_identity(b128mat_t *b) {
  b128mat_t a = B128MAT_INITIALIZER;
  b128_init(b->size,&a);
  bzero(a.b,a.size*a.colb*sizeof(uint128_t));
  madd_identity(&a);
  return a;
}


// matrix addition
// sum size x size b128 compressed matrices :a and :b storing
// the result in c, old content of c is discarded
// :a and :b must be of the same size
// Not tested!
void msum(const b128mat_t *a, const b128mat_t *b, b128mat_t *c)
{
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
void mmult(const b128mat_t *a, const b128mat_t *b, b128mat_t *c)
{
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
// 1) the size of the matrix (a size_t)
// 2) the bit array (an uint128_t array of size a->size*a->colb)
void msave_to_file(const b128mat_t *a, const char *filename)
{
  assert(a!=NULL);
  size_t size = a->size;
  FILE *f = fopen(filename,"w");
  if(f==NULL) quit("msave_to_file: cannot open file", __LINE__,__FILE__);
  size_t e = fwrite(&size,sizeof(size),1,f);
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
size_t mload_from_file(b128mat_t *a, const char *filename)
{
  assert(a!=NULL);
  b128_free(a);
  FILE *f = fopen(filename,"r");
  if(f==NULL) quit("mload_from_file: cannot open file", __LINE__,__FILE__);
  size_t size;
  size_t e = fread(&size,sizeof(size),1,f);
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

// do nothing, added for compatibility with k2 matrices
void minimat_init(int m)
{
  (void) m;
}
void minimat_reset() {}

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
static void b128_free(b128mat_t *m)
{
  if(m->read_only) // read only matrices are pointers to other matrices
    quit("Illegal operation: freeing a read only b128-matrix",__LINE__,__FILE__); 
  if(m->b!=NULL) free(m->b);
  m->b=NULL;
  m->size = m->colb = 0;
}

// make :a a b128 matrix of a given :size by initializing its fields
static void b128_init(size_t size, b128mat_t *a)
{
  assert(a!=NULL);
  if(a->b!=NULL) quit("b128_init: initalizing a non-empty matrix",__LINE__,__FILE__);
  if(size<=0) quit("b128_init: illegal matrix size",__LINE__,__FILE__);
  a->size = size;
  a->colb = (uint32_t) (size+127)/128;
  a->b = malloc(a->size*a->colb*sizeof(uint128_t));
  if(a->b==NULL) quit("b128_init: malloc failed",__LINE__,__FILE__);
  return;
}


