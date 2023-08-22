// functions for bbm matrices
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include "bbm.h"


// write a size x size submatrix containing the value b inside a bbm matrix m
// starting at position i,j entries outsize the matrix m should not be written 
// used for initliazing a matrix or for uncompressing a k2 matrix
void byte_to_bbm(uint8_t *m, int msize, int i, int j, int size, uint8_t b) {
  assert(i>=0 && j>=0 && i<msize+2*size && j<msize+2*size);
  for(int ii=0; ii<size; ii++)
    for(int jj=0; jj<size; jj++)
      if(i+ii<msize && j+jj<msize) 
        m[(i+ii)*msize + j+jj]=b;  // set m[i+ii][j+jj] to b
}


// write the content of a bbm submatrix m to a file f
void bbm_to_ascii(uint8_t *m, int msize, int i, int j, int size, FILE *f)
{
  assert(i>=0 && j>=0 && i<msize && j<msize);
  fprintf(f,"Submatrix at (%d,%d) of size %d\n",i,j,size);  
  for(int ii=0; ii<size; ii++) {
    for(int jj=0; jj<size; jj++) {
      if(i+ii<msize && j+jj<msize) 
        fprintf(f,"%d",m[(i+ii)*msize + j+jj]);  
      else fprintf(f,".");
    }
    fprintf(f,"\n");
  }
}

// multiplication of bbm matrices, return number of nonzero 
// all mustrice must have been allocated to dim size*size
int mmult_bbm(const uint8_t *a, int size, const uint8_t *b, uint8_t *c) {
  assert(a!=NULL && b!=NULL && c!=NULL && size>0);
  int count=0;
  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++) {
      int sum=0;
      for(int k=0; k<size; k++) 
        sum |= a[i*size+k] & b[k*size+j];
      if (sum) {
        c[i*size+j] = 1;
        count++;
      }
    }
  return count;
}

bool mequals_bbm(const uint8_t *a, int size, const uint8_t *b) {
  assert(a!=NULL && b!=NULL && size>0);
  for(int i=0; i<size*size; i++)
      if(a[i] != b[i]) return false;
  return true;
}

