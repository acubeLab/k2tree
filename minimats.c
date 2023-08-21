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
// to initialized in minimats_init();
static minimat_t *mprods2x2 = NULL; // will contain 256 entries  ;

// multiply two minimat matrices of size 2*2
minimat_t mmult2x2(minimat_t a, minimat_t b) {
  assert(mprods2x2!=NULL);
  return mprods2x2[(a<<4) | b]; // equivalent to minimat_prods[a][b]
}

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
      c =  (arow[0] | bcol[0]) ? 1 : 0; // c[0][0] = a[0]*b[0]
      c |= (arow[0] | bcol[1]) ? 2 : 0; // c[0][1] = a[0]*b[1]
      c |= (arow[1] | bcol[0]) ? 4 : 0; // c[1][0] = a[1]*b[0]
      c |= (arow[1] | bcol[1]) ? 8 : 0; // c[1][1] = a[1]*b[1]
      mprods2x2[a<<4 | b] = c;
    }
}

void minimat_init(int msize) {
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
  if(MMsize!=2) quit("minimats_init: MMsize!=2",__LINE__,__FILE__);  
  init_mprods2x2();
}

// compute the size of the smallest k2mat containing a matrix of size msize
// the size depends on the size of the minimat and grows with powers of 2
int k2get_k2size(int msize) 
{
  if(4*Minimat_node_ratio != (MMsize*MMsize)) 
    quit("getk2size: minimats_init not called",__LINE__,__FILE__);
  int s = 2*MMsize;    // size of the smallest legal k2mat
  while(s < msize) {
    s*=2;
    if(s<0) quit("getk2size: overflow",__LINE__,__FILE__);
  }
  if(s>(1<<30)) quit("getk2size: overflow",__LINE__,__FILE__);
  return s;
}

// read a minimat from a submatrix of an compressed matrix m
// entries outsize the matrix m are considered to be 0
minimat_t minimat_read(uint8_t *m, int msize, int i, int j, int size) {
  assert(size==MMsize);
  assert(i>=0 && j>=0 && i<msize+2*size && j<msize+2*size);
  minimat_t res = 0;
  for(int ii=0; ii<size; ii++)
    for(int jj=0; jj<size; jj++) {
      res <<= 1;
      if(i+ii<msize && j+jj<msize && m[(i+ii)*msize + j+jj])
        res |= 1;
    }
  return res;
}



// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
