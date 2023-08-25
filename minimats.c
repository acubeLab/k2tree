/* Routines for arithmetic/logical operations on binary matrices 
   represented as k^2 trees 

   This file contains the functions releated to minimatrices, ie,
   the matrices stored as leaves of the k^2 tree.

   minimat matrices have size MMsize which is set (and never changed)
   at the start of the execution by calling minimats_init
   
   MMsize can be any even value provided a minimat fits a minimat_t 
   variable, and we provide the code for multipliying matrices of size MMsize
   Currently the only supported size is 2 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <string.h>
#include "k2.h"

static void quit(const char *msg, int line, char *file);

// the type representing a (temporary) minimat matrix must be an integer
// so that two minimat matrices can be summed with a single
// bitwise OR operation
typedef uint64_t minimat_t; // explict matrix aka minimat (currently 2x2)

// ----------------------------------------------
// global variable storing the size of a mini-matrix 
// ie. the last level of recursion
// possibly this will be a command line parameter
static int MMsize = INT32_MAX;  // a certainly incorrect value 
// global variable storing how many nodes takes a minimat:
// since nodes have 4 bits (not likely to change)
// this amounts to how many nibbles takes a minimat 
// eg for a 2x2 binary matrix takes one nibble so the constant is 1
static int Minimat_node_ratio = -1;  // certainly incorrect value

// minimat constants (depend on size, here are 2x2)
static minimat_t MINIMAT0s=0;  // minimat containing all 0's, correctly initialized 
static minimat_t MINIMAT1s=0;  // minimat containing all 1's, incorrect value, initialzed in minimats_init  


// multiplication functions for minimats of different sizes

// global variable containing the product of every pair of possible minimatrices
// to be initialized in minimats_init();
static minimat_t *mprods2x2 = NULL; // will contain 256 entries  ;

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
void init_mprods2x2(void) {
  // create array
  assert(mprods2x2==NULL);
  mprods2x2 = malloc(256*sizeof(minimat_t));
  if(mprods2x2==NULL) quit("init_mprods2x2: malloc failed",__LINE__,__FILE__);
  // fill with products
  minimat_t arow[2], bcol[2],c;
  for(minimat_t a=0; a<16; a++) 
    for(minimat_t b=0; b<16;b++) {
      arow[0] = a&3; arow[1] = a>>2;
      bcol[0] = (b&1) | (b&4)>>1; 
      bcol[1] = (b&2)>>1 | (b&8)>>2;
      c =  (arow[0] & bcol[0]) ? 1 : 0; // c[0][0] = a[0]*b[0]
      c |= (arow[0] & bcol[1]) ? 2 : 0; // c[0][1] = a[0]*b[1]
      c |= (arow[1] & bcol[0]) ? 4 : 0; // c[1][0] = a[1]*b[0]
      c |= (arow[1] & bcol[1]) ? 8 : 0; // c[1][1] = a[1]*b[1]
      mprods2x2[a<<4 | b] = c;
    }
}

// init array of transpose
void init_mtranspose4x4(void) {
  assert(mtranspose4x4==NULL);
  mtranspose4x4 = malloc((1<<16)*sizeof(*mtranspose4x4));
  if(mtranspose4x4==NULL) quit("init_mtranspose4x4: malloc failed",__LINE__,__FILE__);
  // fill with transpose
  minimat_t bt[4];      // rows of b^t   
  for(minimat_t b=0; b<(1<<16); b++) {
    for(int j=0;j<4;j++) { // column index in b
      bt[j] = 0;           // the j-th column of b determines the j-th row of bt   
      for(int r=0;r<4;r++) // row index in b
        if(b & (1<<(j+r*4)))  // if b[r][j]!=0
          bt[j] |= (1<<r);
    }
    mtranspose4x4[b] = bt[0] | (bt[1]<<4) | (bt[2]<<8) | (bt[3]<<12);
  }
}



void minimat_init(int msize) {
  if(MMsize!=INT32_MAX) quit("minimats_init: already initialized",__LINE__,__FILE__);
  // msize must be even to ensure that a mimimats has a bit size multiple of 4
  // ie a size multiple of the size of a tree node. 
  if(msize%2!=0) 
    quit("The size of minimat must be even (see code)",__LINE__,__FILE__);
  if(msize*msize>8*sizeof(minimat_t))
    quit("Minimat size too large for type minimat_t",__LINE__,__FILE__);   
        
  // init globals MMsize and Minimat_node_ratio
  MMsize = msize;
  Minimat_node_ratio = (MMsize*MMsize)/4;
 // init MINIMAT1s
 assert(MINIMAT0s==0);
 if(msize*msize == 8*sizeof(minimat_t))
   MINIMAT1s = ~0; // minimat takes the whole variable
 else 
   MINIMAT1s = (((minimat_t) 1) << (msize*msize)) -1;

  // so far only size 2 is allowed
  if(MMsize==2)  
    init_mprods2x2();
  else if (MMsize==4) {
    init_mtranspose4x4();
  }
  else quit("minimats_init: MMsize!=2,4",__LINE__,__FILE__); 
}

// compute the size of the smallest k2mat containing a matrix of size msize
// the size depends on the size of the minimat and grows with powers of 2
int k2get_k2size(int msize) 
{
  if(4*Minimat_node_ratio != (MMsize*MMsize)) 
    quit("k2get_size: minimats_init not called",__LINE__,__FILE__);
  int s = 2*MMsize;    // size of the smallest legal k2mat
  while(s < msize) {
    s*=2;
    if(s<0) quit("getk2size: overflow",__LINE__,__FILE__);
  }
  if(s>(1<<30)) quit("getk2size: overflow",__LINE__,__FILE__);
  return s;
}

// read a minimat from a submatrix of an bbm matrix m
// entries outsize the matrix m are considered to be 0
minimat_t minimat_from_bbm(uint8_t *m, int msize, int i, int j, int size) {
  assert(size==MMsize);
  assert(i>=0 && j>=0 && i<msize+2*size && j<msize+2*size);
  minimat_t res = 0;
  for(int ii=size-1; ii>=0; ii--)
    for(int jj=size-1; jj>=0; jj--) {
      res <<= 1;
      if(i+ii<msize && j+jj<msize && m[(i+ii)*msize + j+jj])
        res |= 1;
    }
  return res;
}

// write the content of a minimat to a bbm matrix m
// entries outsize the matrix m should not be written
void minimat_to_bbm(uint8_t *m, int msize, int i, int j, int size, minimat_t a) {
  assert(size==MMsize);
  assert(i>=0 && j>=0 && i<msize+2*size && j<msize+2*size);
  for(int ii=0; ii<size; ii++)
    for(int jj=0; jj<size; jj++)
      if(i+ii<msize && j+jj<msize) {       // inside the matrix
        int bit = a & (1<<(ii*size+jj));   // read bit a[ii][jj]
        m[(i+ii)*msize + j+jj] = bit ? 1 : 0;  // if nonzero set m[i+ii][j+jj] to 1
        // if m initialized with 0s we can use the following line and save some writes
        // if(bit) m[(i+ii)*msize + j+jj]=1;  // if nonzero set m[i+ii][j+jj] to 1
      }
}






// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
