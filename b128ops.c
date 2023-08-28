/* Routines for arithmetic/logical operations on binary matrices 
   represented as k^2 trees 

   This file contains the defintion of the operations that make use of
   the basic operations defined in b128aux.c

   Matrix dimensions are assumed to be power of 2, of size at least
   2*MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "b128.c"
#include "bbm.h"



// write the content of the :size x :size b128 matrix :a to the bbm matrix :m 
// of size msize*msize. It is assumed m was already correctly initialized 
void mwrite_to_bbm(uint8_t *m, int msize, int size, const b128mat_t *a)
{
  assert(size>=msize);
  if(b128is_empty(a)) {  // an empty b128 matrix is all 0s
    byte_to_bbm(m,msize,0,0,msize,0); // fill m with 0s
    return;
  }
  byte_to_bbm(m,msize,0,0,msize,2); // fill m with illegal value 2
  size_t pos = 0;
  mdecode(m,msize,0,0,size,a,&pos);
  assert(pos==b128pos(a)); // check we read all the b128mat_t structure
}


// compress the matrix *m of size msize into the b128mat_t structure *a 
// m should be an integer array of size msize*msize 
// return the size of the b128 matrix (which has the form 2**k*MMsize)
int mread_from_bbm(uint8_t *m, int msize, b128mat_t *a)
{
  assert(a!=NULL && m!=NULL);
  b128_free(a); // free previous content of a
  assert(msize>1);
  if(msize>(1<<30)) quit("mread_from_bbm: matrix too large",__LINE__,__FILE__);
  int asize = b128get_b128size(msize);
  assert(asize>=2*MMsize);
  // read matrix m into a
  mencode(m,msize,0,0,asize,a);
  return asize;
}

// write to :file statistics for a b128 matrix :a with an arbitrary :name as identifier
void mshow_stats(size_t size, int asize, const b128mat_t *a, const char *mname,FILE *file) {
  (void) asize;
  if(size!=a->size) quit("mshow_stata: size mismatch",__LINE__,__FILE__);
  fprintf(stderr,"%s -- matrix size: %zd, column blocks %d\n",mname,size,a->colb);  
}

// main entry point for matrix addition
// sum size x size b128 compressed matrices :a and :b storing
// the result in c, must be called with an initialized and empty :c
// :a and :b must be of size at least 2*MMsize but their content can be 
// arbitrary: all 0's, all 1's, or generic
// at exit:
//    if the result is a zero matrix c is left empty
//    if the result is an all one's matrix c contains a single ALL_ONES node
//    otherwise c is a node + the recursive description of its subtree  
void msum(int size, const b128mat_t *a, const b128mat_t *b, b128mat_t *c)
{
  assert(size>=2*MMsize);
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(b128is_empty(c) && !c->read_only);

  if(b128is_empty(a) && b128is_empty(b)) 
    return;                 // if a==0 && b==0: c=0
  else if(b128is_empty(b))      
    mcopy(size,a,c);        // if b==0: c=a
  else if(b128is_empty(a))    
    mcopy(size,b,c);        // if a==0: c=b
  else { // a and b are both not empty, call msum_rec
    size_t posa=0,posb=0;
    msum_rec(size,a,&posa,b,&posb,c);
    assert(posa==b128pos(a) && posb==b128pos(b)); // a and b were completeley read
    assert(!b128is_empty(c));  // implied by a+b!=0
  }
  return;
}



// main entry point for matrix equality
// check if size x size b128 compressed matrices :a and :b are equal
// if a==b return -d, where d>0 is the number of levels traversed  
// if a!=b return the level>=0 containing the first difference
// (first in the sense of the first level encountered in dfs order)
// Note that if a==b we return the number of visited levels (tree depth), 
// while if a!=b we return the level of the first difference counting from 0 (root)
// these two values differ by one (these can be seen in the returned value)
// :a and :b must be of size at least 2*MMsize but their content can be
// arbitrary: all 0's, all 1's, or generic
// note: here all 0's matrices are considered of depth 1 even if they are empty
int mequals(int size, const b128mat_t *a, const b128mat_t *b)
{
  assert(size>=2*MMsize);
  assert(a!=NULL && b!=NULL);
  if(b128is_empty(a) && b128is_empty(b)) 
    return -1;                 // if a==0 && b==0: a==b and one level traversed
  else if(b128is_empty(b))      
    return 0;                 // if b==0 && a!=0: a!=b and difference at level 0
  else if(b128is_empty(a))    
    return 0;                 // if a==0 && b!=0: a!=b as above
  // a and b are both not zero
  size_t posa=0,posb=0;
  int eq = mequals_rec(size,a,&posa,b,&posb);
  // do extra checks if the matrices are equal
  assert(eq>=0 || (posa==b128pos(a) && posb==b128pos(b) && posa==posb));
  return eq;
}



// free mem, *m still reusable if needed 
void b128_free(b128mat_t *m)
{
  if(m->read_only) // read only matrices are pointers to other matrices
    quit("Illegal operation: freeing a read only b128-matrix",__LINE__,__FILE__); 
  if(m->b!=NULL) free(m->b);
  m->b=NULL;
  m->size = m->colb = 0;
}

void b128_init(int size, const b128mat_t *a)
{
  assert(a!=NULL);
  assert(!a->read_only);
  if(size<=0) quit("b128_init: illegal matrix size",__LINE__,__FILE__);
  a->size = size;
  a->colb = (size+127)/128;
  a->b = realloc(a->b,a->size*a->colb*sizeof(uint128_t));
  if(a->b==NULL) quit("b128_init: malloc failed",__LINE__,__FILE__);
  return;
}

// main entry point for matrix multiplication. 
// multiply size x size b128 compressed matrices :a and :b storing
// the result in c, must be called with an initialized and empty :c
// :a and :b must be of size at least 2*MMsize but their content can be 
// arbitrary: all 0's, all 1's, or generic
// at exit:
//    if the result is a zero matrix c is left empty
//    if the result is an all one's matrix c contains a single ALL_ONES node
//    otherwise c is a node + the recursive description of its subtree  
void mmult(int asize, const b128mat_t *a, const b128mat_t *b, b128mat_t *c)
{
  (void) asize;
  assert(a!=NULL && b!=NULL && c!=NULL);
  if(a->size != b->size)
    quit("mmult: matrix size mismatch",__LINE__,__FILE__);

  // remove????    
  //assert(b128is_empty(c));
  //if one matrix is all 0s (empty) the result is all 0's: nothing to be done  
  //if(b128is_empty(a) ||  b128is_empty(b))
  //  return;

  // init and clean c
  uint128_t one = 1;
  b128_init(a->size,c);
  assert(c->colb==a->colb && c->colb == b->colb);
  bzero(c->b,c->size*c->colb*sizeof(uint128_t));
  for(size_t i=0; i<a->size; i++)  // row index for a 
    for(size_t k=0; k<a->size; k++) { // col index for a 
      size_t pos = i*a->colb+k/128;  
      if(a->b[pos] & (one << (k%128) )) { // a[i][k]==1
        // c[i,:] != b[k,:]
        for(int j=0; j<a->colb; j++)  // iterate over all row blocks 
          c[i*a->colb + j] |= b[k*a->colb +j];
      }
    }
}

// return statistics on matrix a
// write 0 in the variables passed by reference
// and return the number of block columns
int mstats(int asize, const b128mat_t *a, size_t *pos, size_t *nodes, size_t *minimats)
{
  void (asize);
  assert(a!=NULL);
  *pos=*nodes=*minimats=0;
  return a->colb;
}


// save the matrix :a to file :filename
// the format is
// 1) the actual size of the matrix (an int)
// 5) the bit array (an uint128_t array of size a->size*a->colb)
void msave_to_file(int size, int asize, const b128mat_t *a, const char *filename)
{
  (void) asize;
  assert(a!=NULL);
  if(size!=a->size) quit("msave_to_file: matrix size mismatch", __LINE__,__FILE__);
  FILE *f = fopen(filename,"w");
  if(f==NULL) quit("msave_to_file: cannot open file", __LINE__,__FILE__);
  size_t e = fwrite(&size,sizeof(int),1,f);
  if(e!=1) quit("msave_to_file: cannot write size",__LINE__,__FILE__);
  size_t tot128 = ((size_t) size)*a->colb;
  if(tot128 <=0) quit("msave_to_file: illegal bit array",__LINE__,__FILE__);
  e = fwrite(&a->b,sizeof(uint128_t),tot128,f);
  if(e!=tot128) quit("msave_to_file: cannot write bit rows",__LINE__,__FILE__);
  fclose(f);
}


// load a b128 matrix stored in file :filename into the b128mat_t structure :a
// return the actual size of the matrix, see msave_to_file() for the file format
// must be called after minimat_init() since it checks that the correct MMsize is used
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
  if(size<1) quit("mload_from_file: matrix size smaller than 2, wrong format?",__LINE__,__FILE__);
  b128_init(size,a);
  size_t = tot128 = a->size*a->colb;
  e = fread(a->b,sizeof(uint128_t),tot128,f);
  if(e!=tot128) quit("mload_from_file: cannot read positions",__LINE__,__FILE__);
  fclose(f);
  return size;
}




// ----------- static auxiliary functions ------------

