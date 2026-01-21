/* Routines for arithmetic/logical operations on binary matrices 
   represented as k^2 trees 

   This file contains the functions related to minimatrices, ie,
   the matrices stored as leaves of the k^2 tree.

   minimat matrices have size MMsize which is set (and never changed)
   at the start of the execution by calling minimats_init
   
   MMsize can be any even value provided a minimat fits a minimat_t 
   variable, and we provide the code for multipliying matrices of size MMsize
   Currently the only supported size is 2 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _K2MINIMATS_C
#define _K2MINIMATS_C 1

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <string.h>
#include "k2.h"
#include "bbm.h"

static void quit(const char *msg, int line, char *file);

// the type representing a (temporary) minimat matrix must be an integer
// so that two minimat matrices can be summed with a single
// bitwise OR operation
typedef uint64_t minimat_t; // explict matrix aka minimat (currently 2x2)

// ----------------------------------------------
// global variable storing the size of a mini-matrix 
// ie. the last level of recursion
// possibly this will be a command line parameter
static uint32_t MMsize = 0;  // a certainly incorrect value 
// global variable storing how many nodes takes a minimat:
// since nodes have 4 bits (not likely to change)
// this amounts to how many nibbles takes a minimat 
// eg for a 2x2 binary matrix takes one nibble so the constant is 1
uint32_t Minimat_node_ratio = 0;  // certainly incorrect value

// global variable determining whether during construction we use the 
// ALL_ONES node to denote a submatrix of all 1's.
bool Use_all_ones_node = true;

// minimat constants (depend on size, these are for 2x2 and 4x4)
static minimat_t MINIMAT0s=0;  // minimat containing all 0's, correctly initialized 
static minimat_t MINIMAT1s;    // minimat containing all 1's, initialized in minimats_init  
static minimat_t MINIMAT_Id;   // identity minimat initialized in minimats_init  

// multiplication functions for minimats of different sizes

// global variable containing the product of every pair of possible minimatrices
// to be initialized in minimats_init();
static minimat_t *mprods2x2 = NULL; // will contain 256 entries
static uint16_t mtranspose2x2[16]; // transpose of each 2x2 minimat
static uint16_t *mtranspose=NULL;  // transpose of minimats, must be initialized in minimats_init

// multiply two minimat matrices of size 2*2
minimat_t mmult2x2(minimat_t a, minimat_t b) {
  assert(mprods2x2!=NULL);
  return mprods2x2[(a<<4) | b]; // equivalent to minimat_prods[a][b]
}

// global variable containing the transpose of each 4x4 minimat 
static uint16_t *mtranspose4x4 = NULL;

minimat_t mmult4x4(minimat_t a, minimat_t b) {
  // compute array of b columns
  minimat_t bt = mtranspose4x4[b];
  minimat_t bcol[4];
  for(int j=0;j<4;j++) {
    bcol[j] = bt & 0xF;
    bt = bt >> 4;
  }
  // compute result matrix
  minimat_t res = 0;
  for(int i=3;i>=0;i--) { // row index
    minimat_t rowi = (a>>(4*i)) & 0xF;
    for(int j=3; j>=0;j--) {
      res <<= 1;
      if( rowi & bcol[j] ) res |= 1;
    }
  }
  return res;
}

// init array of products
static void init_mprods2x2(void) {
  // create array
  assert(mprods2x2==NULL);
  mprods2x2 = malloc(256*sizeof(minimat_t));
  if(mprods2x2==NULL) quit("init_mprods2x2: malloc failed",__LINE__,__FILE__);
  // fill with products
  minimat_t arow[2], bcol[2],c;
  for(minimat_t a=0; a<16; a++) {
    for(minimat_t b=0; b<16;b++) { // init products times b
      arow[0] = a&3; arow[1] = a>>2;
      bcol[0] = (b&1) | (b&4)>>1; 
      bcol[1] = (b&2)>>1 | (b&8)>>2;
      c =  (arow[0] & bcol[0]) ? 1 : 0; // c[0][0] = a[0]*b[0]
      c |= (arow[0] & bcol[1]) ? 2 : 0; // c[0][1] = a[0]*b[1]
      c |= (arow[1] & bcol[0]) ? 4 : 0; // c[1][0] = a[1]*b[0]
      c |= (arow[1] & bcol[1]) ? 8 : 0; // c[1][1] = a[1]*b[1]
      mprods2x2[a<<4 | b] = c;
    }
    // init transpose array
    mtranspose2x2[a] = (uint16_t) ( (a&1) | ((a&2)<<1) | ((a&4)>>1) | (a&8) );
  }
}

// init array of transpose
static void init_mtranspose4x4(void) {
  assert(mtranspose4x4==NULL);
  assert(MMsize==4);
  mtranspose4x4 = malloc((1<<16)*sizeof(*mtranspose4x4));
  if(mtranspose4x4==NULL) quit("init_mtranspose4x4: malloc failed",__LINE__,__FILE__);
  // fill with transpose
  minimat_t bt[4];      // rows of b^t   
  for(minimat_t b=0; b<(1<<16); b++) {
    for(int j=0;j<4;j++) { // column index in b
      bt[j] = 0;           // the j-th column of b determines the j-th row of bt   
      for(int r=0;r<4;r++) // row index in b
        if(b & (1<<(j+r*4)))  // if b[r][j]!=0
          bt[j] |= (1UL<<r);
    }
    mtranspose4x4[b] = (uint16_t) (bt[0] | (bt[1]<<4) | (bt[2]<<8) | (bt[3]<<12));
  }
}


// init minimat constants
void minimat_init(int msize) {
  if(MMsize!=0) quit("minimats_init: already initialized",__LINE__,__FILE__);
  // msize must be even to ensure that a mimimats has a bit size multiple of 4
  // ie a size multiple of the size of a tree node. 
  if(msize%2!=0) 
    quit("The size of minimat must be even (see code)",__LINE__,__FILE__);
  if(msize*msize>8*sizeof(minimat_t))
    quit("Minimat size too large for type minimat_t",__LINE__,__FILE__);   
        
  // init globals MMsize and Minimat_node_ratio
  MMsize = (uint32_t) msize;
  Minimat_node_ratio = (MMsize*MMsize)/4;
 // init MINIMAT1s
 assert(MINIMAT0s==0);
 if(msize*msize == 8*sizeof(minimat_t))
   MINIMAT1s = (minimat_t) ~0; // minimat takes the whole variable
 else 
   MINIMAT1s = (((minimat_t) 1) << (msize*msize)) -1;

  // so far only size 2 or 4 are allowed
  if(MMsize==2) {
    init_mprods2x2();
    mtranspose = mtranspose2x2;
    MINIMAT_Id = 0xA; // TODO check this!
  }
  else if (MMsize==4) {
    init_mtranspose4x4();
    mtranspose = mtranspose4x4;
    MINIMAT_Id = 0x8421;
  }
  else quit("minimats_init: MMsize!=2,4",__LINE__,__FILE__); 
}

// return minimat size
uint32_t minimat_size() {
  if(MMsize==0) quit("minimats_size: not initialized",__LINE__,__FILE__);
  return MMsize;
}


// reset minimats initialization and makes it possible
// to call minimat init again
void minimat_reset() {
  MMsize=0;
  Minimat_node_ratio=0;
  mtranspose = NULL;
  if(mprods2x2!=NULL) {
    free(mprods2x2);
    mprods2x2=NULL;
  }
  if(mtranspose4x4!=NULL) {
    free(mtranspose4x4);
    mtranspose4x4=NULL;
  }
}


// read a minimat from a submatrix of a bbm matrix m of size msize
// entries outsize the matrix m are considered to be 0
minimat_t minimat_from_bbm(uint8_t *m, size_t msize, size_t i, size_t j, size_t size) {
  // printf("i=%zu j=%zu size=%zu\n",i,j,size);  
  assert(size==MMsize);
  assert(i<msize+2*size && j<msize+2*size);
  minimat_t res = 0;
  for(ssize_t ii=size-1; ii>=0; ii--)
    for(ssize_t jj=size-1; jj>=0; jj--) {
      res <<= 1;
      if(i+ii<msize && j+jj<msize && m[(i+ii)*msize + j+jj])
        res |= 1;
    }
  return res;
}

// read a minimat from the interleaved array ia[0,n-1] containing entries in [imin,imax)
minimat_t minimat_from_ia(uint64_t *ia, size_t n, uint64_t imin, size_t size) {
  assert(size==MMsize);     // only called for minimats
  assert(n>0);              // not called on an empty submatrix  
  assert(n<=MMsize*MMsize); // cannot have more than MMsize*MMsize entries
  minimat_t res = 0;
  if(MMsize==2) {
    // assert(imax==imin+4);    imax removed from parameters
    for(size_t i=0; i<n; i++) {
      assert(ia[i]>=imin && ia[i]-imin <4); // equiv to imin <= ia[i]<imax
      int64_t j = (int64_t) (ia[i]-imin); 
      // for MMsize=2 j is the position of the corresponding 1 in res
      assert(j>=0 && j<4);
      res |= (1UL<<j);
    }
  }
  else if(MMsize==4) {
    // map between entries in interleaved vs row-major order for size==4
    int t[16] = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};
    // assert(imax==imin+16);    imax removed from parameters
    for(size_t i=0; i<n; i++) {
      assert(ia[i]>=imin && ia[i]-imin <16); // equiv to imin <= ia[i]<imax
      int j = t[ia[i]-imin]; 
      assert(j>=0 && j<16);
      res |= (1UL<<j);
    }
  }
  else quit("minimat_from_ia: MMsize!=2,4",__LINE__,__FILE__); 
  assert(res!=MINIMAT0s); // cannot be all 0's
  return res;  
}


// write the content of a minimat of size size 
// to a bbm matrix m of size msize
// entries outsize the matrix m should not be written
void minimat_to_bbm(uint8_t *m, size_t msize, size_t i, size_t j, size_t size, minimat_t a) {
  assert(size==MMsize);
  assert(i<msize+2*size && j<msize+2*size);
  for(size_t ii=0; ii<size; ii++)
    for(size_t jj=0; jj<size; jj++)
      if(i+ii<msize && j+jj<msize) {       // inside the matrix
        minimat_t bit = a & (1<<(ii*size+jj));   // read bit a[ii][jj]
        m[(i+ii)*msize + j+jj] = bit ? 1 : 0;  // if nonzero set m[i+ii][j+jj] to 1
        // if m initialized with 0s we can use the following line and save some writes
        // if(bit) m[(i+ii)*msize + j+jj]=1;  // if nonzero set m[i+ii][j+jj] to 1
      }
}

// write the content of a minimat to a text file f
// as a list of arcs. there should be no arcs outside msizeXmsize
// msize is the size of the orginal matrix
// i,j is the upper left corner of the minimat
// :size is the size of the minimat :a
void minimat_to_text(FILE *f, size_t msize, size_t i, size_t j, size_t size, minimat_t a) {
  assert(size==MMsize);
  assert(i<msize+2*size && j<msize+2*size);
  for(size_t ii=0; ii<size; ii++)
    for(size_t jj=0; jj<size; jj++)
      if(i+ii<msize && j+jj<msize) {       // inside the matrix
        minimat_t bit = a & (1<<(ii*size+jj));   // read bit a[ii][jj]
        if(bit) { 
          int e = fprintf(f,"%zu %zu\n",i+ii,j+jj);
          if(e<0) quit("minimat_to_text: fprintf failed",__LINE__,__FILE__);
        }
      }
      else assert((i==j && ii==jj) ||  ((a & (1<<(ii*size+jj))) ==0) ); // no entry outside msize except the diagonal
}


// matrix-vector multiplication for a minimat matrix of size 2x2
// inside the minimat the bit matrix is stored in row-major order
// a_00 is the less significant bit, a_01 the next one, then a_10 and a_11
void mvmult2x2(minimat_t a, const vfloat *x, vfloat *y) {
  assert(a!=MINIMAT0s); // a must be non zero (tested before calling the function)
  assert(a!=0);         // same condition as above but more explicit
  if (a & 1) y[0] += x[0];
  if (a & 2) y[0] += x[1];
  if (a & 4) y[1] += x[0];
  if (a & 8) y[1] += x[1];
  return;
}

// matrix-vector multiplication for a minimat matrix of size 4x4
// inside the minimat the bit matrix is stored in row-major order
// a_00 is the less significant bit, then a_01 a_02 a_03 a_10, a_11 ...
void mvmult4x4(minimat_t a, const vfloat *x, vfloat *y) 
{
  assert(a!=MINIMAT0s); // a must be non zero (tested before calling the function)
  assert(a!=0);         // same condition as above but more explicit
  for(size_t i=0; i<4 && a!=0; i++) {
    minimat_t row = a & 0xF;
    if(row & 1) y[i] += x[0];
    if(row & 2) y[i] += x[1];
    if(row & 4) y[i] += x[2];
    if(row & 8) y[i] += x[3];
    a >>= 4;
  }
  return;
}



// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
#endif