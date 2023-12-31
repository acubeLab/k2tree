/* functions for bbm matrices (binary byte matrix)

   Technical note: matrix sizes are int, therefore limited to 2^31, but values related
   to the overall number of elements can be up to 2^61 therefore are always 
   stored into a size_t variable (usually 64 bits))  

   Hence a size x size bbm matrix is represented by a one-dimensional uint8_t array b[]
   and of its length, equals to size*size, stored in a size_t variable
   The entries should be 0/1, ie two different nonzero values are considered 
   different by mequals_bbm(). 
*/
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include "bbm.h"


// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

// compute integer square root
static size_t intsqrt(size_t n) {
  // assert(n>=0);
  size_t x = n;
  size_t y = (x + 1) / 2;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2;
  }
  assert(x*x <= n && (x+1)*(x+1) > n);
  return x;
}

// alloc memory for a size x size bbm matrix
uint8_t *bbm_alloc(size_t size)
{
  size_t length = size*size;
  uint8_t *buffer = malloc(length);
  if(buffer==NULL) quit("Out of memory",__LINE__,__FILE__);
  return buffer;
}

void bbm_write(const uint8_t *m, size_t msize, const char *name)
{
  FILE *f = fopen(name,"wb");
  if(f==NULL) 
    quit("Cannot open matrix file",__LINE__,__FILE__);
  size_t w = fwrite(m, 1, msize*msize, f);
  if(w!=msize*msize) 
    quit("Cannot write matrix file",__LINE__,__FILE__);
  fclose(f);
}


// given a file name of a matrix in bbm format
// return byte array and its size
// exit program on error 
uint8_t *bbm_read(const char *name, size_t *psize)
{
  FILE *f = fopen(name,"rb");
  if(f==NULL) 
    quit("Cannot open matrix file",__LINE__,__FILE__);
  // get file size
  fseek(f, 0, SEEK_END);
  size_t length = ftell(f);
  // check if square
  size_t size = intsqrt(length);
  if(size*size != length) 
    quit("Non square input matrix",__LINE__,__FILE__);  

  // save matrix size
  *psize = size;
  // read file into buffer
  rewind(f);
  uint8_t *buffer = malloc(length);
  if(buffer==NULL) quit("Out of memory",__LINE__,__FILE__);
  size_t r = fread(buffer, 1, length, f);
  if(r!=length) quit("Cannot read matrix file",__LINE__,__FILE__);  
  fclose(f);
  return buffer;
}



// write a size x size submatrix containing the value b inside a bbm matrix m
// starting at position i,j entries outsize the matrix m should not be written 
// used for initliazing a matrix or for uncompressing a k2 matrix
void byte_to_bbm(uint8_t *m, size_t msize, size_t i, size_t j, size_t size, uint8_t b) {
  assert(i<msize+2*size && j<msize+2*size);
  for(size_t ii=0; ii<size; ii++)
    for(size_t jj=0; jj<size; jj++)
      if(i+ii<msize && j+jj<msize) 
        m[(i+ii)*msize + j+jj]=b;  // set m[i+ii][j+jj] to b
}


// write the content of a bbm submatrix m to a file f
void bbm_to_ascii(const uint8_t *m, size_t msize, size_t i, size_t j, size_t size, FILE *f)
{
  assert(i>=0 && j>=0 && i<msize && j<msize);
  fprintf(f,"Submatrix at (%d,%d) of size %d\n",i,j,size);  
  for(size_t ii=0; ii<size; ii++) {
    for(size_t jj=0; jj<size; jj++) {
      if(i+ii<msize && j+jj<msize) 
        fprintf(f,"%d",m[(i+ii)*msize + j+jj]);  
      else fprintf(f,".");
    }
    fprintf(f,"\n");
  }
}

// multiplication of bbm matrices, return 0 (it used to return the 
// number of nonzero elements in the product but this is no longer available) 
// all matrices must have been allocated to dim size*size
void mmult_bbm(const uint8_t *a, size_t size, const uint8_t *b, uint8_t *c) {
  assert(a!=NULL && b!=NULL && c!=NULL && size>0);
  // clean c
  byte_to_bbm(c,size,0,0,size,0);
  // different access pattern 
  for(size_t i=0; i<size; i++)
    for(size_t k=0; k<size; k++) 
      if(a[i*size+k]) {
        for(size_t j=0; j<size; j++) { 
          if(b[k*size+j])
            c[i*size+j] = 1;
        }
      }

  #if 0
  // this is very very slow because of its memory access pattern
  int count=0; // remove this variable in future versions
  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++) {
      // exit for loop as soon as we know result is one
      for(int k=0; k<size; k++) 
        if( a[i*size+k] & b[k*size+j] ) {
          c[i*size+j] = 1;
          count++;
          break;
        }
    }
  return count;
  #endif
}

// compare two bbm matrices of the same size, return the position 
// of the first difference or -1 if they are equal 
ssize_t mequals_bbm(const uint8_t *a, size_t size, const uint8_t *b) {
  assert(a!=NULL && b!=NULL && size>0);
  for(size_t i=0; i<size*size; i++)
      if(a[i] != b[i]) return i;
  return -1;
}

