/* Routines for arithmetic operations on binary matrices represented as k^2 trees 
   operations are defined with the or/and operators playing the role of sum/product

   This file contains the definitions of the complex operations that make 
   use of the basic operations defined in k2aux.c

   Matrix dimensions are assumed to be power of 2, of size at least
   2*MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#include "rank_0000.h"
#endif
#pragma GCC target ("sse4.2")  // ensure sse4.2 compiler switch it is used 
#include "k2aux.c"    // includes minimats.c k2.h bbm.h
#include "k2text.c"   // functions for input/output of k2 matrices from sparse text files  
#include "k2io.c"     // input/output of of k2 files in compressed format

/* The operations defined in this file are:
  
    mequals: check if two k2 matrices are equal
    msum: sum of two k2 matrices
    mmult: multiplication of two k2 matrices
    mvmult: multiplication of a k2 matrix by a vector
  
    All operations assume that the input matrices are of size at least
    2*MMsize (minimatrix size), ie, the size of the last level of recursion.
  
    The operations are defined recursively, splitting the input matrices
    into four quadrants at each level, until reaching the minimat level.
    The basic operations on nodes and minimats are defined in k2aux.c 
    and minimats.c respectively.
  
    The zero matrix is represented by an empty tree (k2is_empty) 

    The results are always in pdf (plain depth first) format, ie, no subtree
    size information is computed: if necessary it can be added with k2dfs_sizes()
    Backpointers also are not computed, but they cannot be easily added.
    
    If the flag Use_all_ones is set, then all ones submatrix compression
    is used in the output matrix. 


    Proposal: add a use_backp flag if it is not set, then 0000 is an all 1s submatrix
    Proposal: add a transpose flag and a main_diag flag.


*/
static void split_and_rec(size_t size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c);
static void mvmult_rec(size_t size, const k2mat_t *a, vfloat *x, vfloat *y);
static void mdecode_and_multiply(size_t size, const k2mat_t *c, size_t *pos, vfloat *x, vfloat *y);
// global variable to force computation of subtree info on the fly
// chenge this to a per-matrix property
bool Extended_edf = false; // compute subtree info on the fly


// recursive test for equality of two k2 matrices both nonzero
// see mequals for details
static int mequals_rec(size_t size, const k2mat_t *a, size_t *posa, 
                          const k2mat_t *b, size_t *posb)
{
  assert(size>=2*MMsize);
  assert(a!=NULL && b!=NULL);
  assert(!k2is_empty(a) && !k2is_empty(b));
  // read root nodes of a and b
  node_t roota = k2read_node(a,*posa); *posa +=1;
  node_t rootb = k2read_node(b,*posb); *posb +=1;
  if(roota!=rootb) return 0; // if root nodes are different: a!=b at this level
  if(roota==ALL_ONES)        // implies rootb == ALL_ONES 
    return -1; // if root nodes are both all 1's: a==b and there are no other levels 
  // root nodes are equal and have children, check children recursively
  if(size==2*MMsize) { // children are minimat matrices
    minimat_t ax[2][2], bx[2][2];
    k2split_minimats(a,posa,roota,ax);
    k2split_minimats(b,posb,rootb,bx);
    for(int k=0;k<4;k++) 
      if(ax[k/2][k%2]!=bx[k/2][k%2]) return 1; // if corresponding minimats are different: a!=b
    return -2; // all minimats are equal: a==b, traversed 2 levels 
  }
  else { // size>2*MMsize: children are k2 matrices, possibly use recursion
    int eq = 0;
    for(int k=0;k<4;k++) {
      if (roota & (1 << k)) {
        assert(rootb & (1 << k)); 
        int eqr = mequals_rec(size/2,a,posa,b,posb);
        if(eqr>0) return eqr+1; // a and b are different at level eqr+1
        if(eqr < eq) eq = eqr;  // keep track of deepest level reached
      }
    }
    return eq -1; // increase depth by one 
  }
}


// main entry point for matrix equality
// check if size x size k2 compressed matrices :a and :b are equal
// if a==b return -d, where d>0 is the number of levels traversed  
// if a!=b return the level>=0 containing the first difference
// (first in the sense of the first level encountered in dfs order)
// Note that if a==b we return the number of visited levels (tree depth), 
// while if a!=b we return the level of the first difference counting from 0 (root)
// these two values differ by one (this can be seen in the returned value)
// :a and :b must be of size at least 2*MMsize but their content can be
// arbitrary: all 0's, all 1's, or generic
// note: here all 0's matrices are considered of depth 1 even if they are empty
int mequals(size_t size, const k2mat_t *a, const k2mat_t *b)
{
  assert(size>=2*MMsize);
  assert(a!=NULL && b!=NULL);
  if(k2is_empty(a) && k2is_empty(b)) 
    return -1;                 // if a==0 && b==0: a==b and one level traversed
  else if(k2is_empty(b))      
    return 0;                 // if b==0 && a!=0: a!=b and difference at level 0
  else if(k2is_empty(a))    
    return 0;                 // if a==0 && b!=0: a!=b as above
  // a and b are both not zero
  size_t posa=0,posb=0;
  int eq = mequals_rec(size,a,&posa,b,&posb);
  // do extra checks if the matrices are equal
  assert(eq>=0 || (posa==k2pos(a) && posb==k2pos(b) && posa==posb));
  return eq;
}

// full copy of matrix a: to b: used instead of sum when one of the matrices is all 0s
//??? TODO: Instead of using recursion, can we just copy the content of a (from offset to pos) to b?
static void mcopy_full(size_t size, const k2mat_t *a, k2mat_t *b)
{
  assert(size>MMsize);  // required by k2copy_rec 
  assert(a!=NULL && b!=NULL);
  assert(!k2is_empty(a));
  assert(!b->read_only);
  assert(!a->read_only);
  assert(!a->open_ended);
  // new version in which the complete content of a from a->offset to a->pos is copied to b
  // ok since the size of the minimats is an integral multiple of the size of a node  
  for(size_t pos=a->offset;pos < a->pos; pos++) {
    node_t n = k2read_node(a,pos);
    k2add_node(b,n);
  }
  // old version based on recursion, could work if a is open-ended
  //size_t pos = 0;
  // k2copy_rec(size,a,&pos,b);
 }

// recursive sum of two k2 (sub)matrices
// assume plain matrices: no subtree info, no backpointers
// a and b must not be all 0s (sub)matrices
//   (if a or b is all 0s this function is not called because the sum is a copy) 
// a and b can be all 1's 
// the output matrix c is normalized as ususal
//  if c is all 0s nothing is written (this should not happen because the sum is an OR)
//  if c is all 1s just the root ALL_ONES is written
//  otherwise we follow the standard rule: root node + subtrees in DFS order
// If Use_all_ones_node is false and a and b do not contain ALL_ONES nodes
// then neither c does   
static void msum_rec_plain(size_t size, const k2mat_t *a, size_t *posa, 
                         const k2mat_t *b, size_t *posb, k2mat_t *c)
{
  assert(size>=2*MMsize); // there must be the root node 
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(!k2submatrix_empty(a,*posa) && !k2submatrix_empty(b,*posb));
  // take care of all 1s matrices: read root without advancing positions
  node_t roota = k2read_node(a,*posa);
  node_t rootb = k2read_node(b,*posb);
  if(roota==ALL_ONES) {
    k2add_node(c,ALL_ONES); 
    *posa+=1; // reach the end of a
    k2dfs_visit_fast(size,b,posb); //scan but ignore b content (subtree info could speedup?)
    return;
  }
  else if(rootb==ALL_ONES) {
    k2add_node(c,ALL_ONES); 
    *posb+=1; // same as above with a and b swapped
    k2dfs_visit_fast(size,a,posa); //scan but ignore a content
    return;
  }
  assert(roota!=ALL_ONES && rootb!=ALL_ONES);
  // a and b are not all 1s, merge children
  *posa += 1; *posb += 1; // skip root nodes already read
  node_t rootc = roota | rootb; // root node of c, correct except when c is all 1s
  assert(rootc!=NO_CHILDREN);   // at least one child is nonzero
  size_t rootpos = k2add_node(c,rootc);  // save c root and its position
  bool all_ones=true;     // true if all submatrices cx[i][j] are all 1's
  if(size==2*MMsize) {    // children are minimat matrices
    minimat_t ax[2][2], bx[2][2];
    k2split_minimats(a,posa,roota,ax);
    k2split_minimats(b,posb,rootb,bx);
    for(int k=0;k<4;k++) {
      minimat_t cx = ax[k/2][k%2] | bx[k/2][k%2]; // compute bitwise or of corresponding minimat
      if (cx != MINIMAT0s) { // save cx if nonzero
        assert(rootc & (1 << k)); // implied by  cx = ax | bx 
        k2add_minimat(c, cx);
      }
      else assert((rootc & (1 << k)) == 0);
      if (cx != MINIMAT1s) all_ones = false;
    }
    // redundant check: all_ones=> 4 minimats stored
    assert(!all_ones || (rootc==ALL_CHILDREN&&k2pos(c)==rootpos+1+4*Minimat_node_ratio)); 
  }
  else { // size>2*MMsize: children are k2 matrices, possibly use recursion
    // we could split a and b in 4 submatrices and sum them, but it is more efficient 
    // to compute the sum without building submatrices (which requires a scan of a and b)
    for(int k=0;k<4;k++) {
      size_t tmp = k2pos(c);        // save current position
      if (roota & (1 << k)) {
        if(rootb & (1 << k)) 
          msum_rec_plain(size/2,a,posa,b,posb,c); // k-th child of c is sum of kth children of a and b
        else 
          k2copy_rec(size/2,a,posa,c); // k-th child of c is kth child of a
      }
      else if (rootb & (1 << k))
        k2copy_rec(size/2,b,posb,c); // k-th child of c is kth child of b
      else 
        assert( (rootc & (1 << k)) == 0); // both children are 0s, nothing to do
      // check tmp and update all_ones
      if(tmp!=k2pos(c)) { // something was written
        assert(k2pos(c)>tmp);
        if(k2read_node(c,tmp)!=ALL_ONES) all_ones = false;
        else assert(k2pos(c)==tmp+1); // the written submatrix was ALL_ONES
      }
      else all_ones = false; // nothing was written, submatrix is all 0s, all_ones is false
    } // end for k=0..3
    assert(!all_ones  || (rootc==ALL_CHILDREN && k2pos(c)==rootpos+5));
  }
  // normalize if c is all 1s (regardless of size)
  if(all_ones && Use_all_ones_node) {
    assert(rootc==ALL_CHILDREN);
    k2setpos(c,rootpos+1);       // discard current children
    k2write_node(c,rootpos,ALL_ONES); // write ALL_ONES as root
  }
  return; 
}

// main entry point for matrix addition
// assume plain matrices: no subtreeinfo no backpointers
// sum size x size k2 compressed matrices :a and :b storing
// the result in :c, the old content of :c is discarded
// :a and :b must be of (same) size, at least 2*MMsize, but their content can be 
// arbitrary: all 0's, all 1's, or generic
// at exit:
//    if the result is a zero matrix, c is left empty
//    if the result is an all one's matrix, c contains a single ALL_ONES node
//    otherwise c is a node + the recursive description of its subtree (note Use_all_ones is ignored) 
// Note: this function is called by matrix product, to sum partial products 
// so the operands are never compressed  
void msum(size_t size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c)
{
  assert(size>=2*MMsize);
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(a->backp==NULL && a->subtinfo==NULL);
  assert(b->backp==NULL && b->subtinfo==NULL);

  k2_free(c); // free old content and initialize as empty
  if(k2is_empty(a) && k2is_empty(b)) 
    return;                 // if a==0 && b==0: c=0
  else if(k2is_empty(b))      
    mcopy_full(size,a,c);        // if b==0: c=a
  else if(k2is_empty(a))    
    mcopy_full(size,b,c);        // if a==0: c=b
  else { // a and b are both not empty, call msum_rec
    size_t posa=0,posb=0;
    msum_rec_plain(size,a,&posa,b,&posb,c);
    assert(posa==k2pos(a) && posb==k2pos(b)); // a and b were completeley read
    assert(!k2is_empty(c));  // implied by a+b!=0
  }
  return;
}

// base case of matrix multiplication: matrices of size 2*MMmin
// a and b must be not all 0s (ie empty)
//   (if a or b is 0s this function is not called because the product is 0) 
// a and b can be all 1's 
// the output matrix c is normalized as usual:
//  if c is all 0s nothing is written
//  if c is all 1s the root ALL_ONES is written
//  otherwise we follow the standard rule: root node + nonzero minisize matrices 
// Note, even if a or b are compressed, there cannot be pointers at this level (height=1)
// Here is the only part where we call the base multiplication function
// using the following macro, change it to support additional sizes
#define mmultNxN(s,a,b) ((s)==2 ? mmult2x2((a),(b)) : mmult4x4((a),(b)))
static void mmult_base(size_t size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c)
{
  assert(size==2*MMsize);
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(!k2is_empty(a) && !k2is_empty(b));
  assert(!c->read_only);
  // c is always an empty matrix because a partial product is never 
  // written directly to a result matrix 
  assert(k2is_empty(c));
  // initialize ax[][] and bx[][] to cover the case when the matrix a/b is all 1s                       
  minimat_t ax[2][2] = {{MINIMAT1s,MINIMAT1s},{MINIMAT1s,MINIMAT1s}};
  minimat_t bx[2][2] = {{MINIMAT1s,MINIMAT1s},{MINIMAT1s,MINIMAT1s}};
  node_t roota = k2read_node(a,0);
  node_t rootb = k2read_node(b,0);
  //both matrices are all 1s ?
  if(roota==ALL_ONES && rootb == ALL_ONES) {
    #ifndef NDEBUG
    if(a->backp!=NULL)
      quit("Illegal left operand: compressed and with an ALL_ONES node",__LINE__,__FILE__);
    if(b->backp!=NULL)
      quit("Illegal right operand: compressed and with an ALL_ONES node",__LINE__,__FILE__);
    #endif
    if(Use_all_ones_node) k2add_node(c,ALL_ONES); // remove this 
    else {
      // if Use_all_ones_node is false we write a 2x2 matrix of all 1s
      k2add_node(c,ALL_CHILDREN); // write root node
      k2add_minimat(c,MINIMAT1s); // write 4 submatrices matrix of all 1s
      k2add_minimat(c,MINIMAT1s);
      k2add_minimat(c,MINIMAT1s);
      k2add_minimat(c,MINIMAT1s);
    }
    return;
  }
  // TODO: insert possible code for case when one matrix is all 1s and the other is not
  // split a and b
  size_t posa=1,posb=1; // we have already read the root node
  if(roota!=ALL_ONES)   // case ALL_ONES is covered by initialization above
    k2split_minimats(a,&posa,roota,ax);
  else assert(a->backp==NULL); // if a is ALL_ONES it cannot have backp pointers 
  if(rootb!=ALL_ONES) 
    k2split_minimats(b,&posb,rootb,bx);
  else assert(b->backp==NULL); // if b is ALL_ONES it cannot have backp pointers

  // split done, now multiply and store c 
  // optimization on all 1's submatrices still missing 
  bool all_ones=true; // true if all c submatrices cx[i][j] are all 1's
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_NODES as root placeholder 
  node_t rootc=NO_CHILDREN;        // actual root node to be computed
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(int k=0;k<4;k++) {  
    int i=k/2; int j=k%2;
    assert(size==4 || size==8); // implies that minimats are 2x2  or 4x4
    minimat_t cx  = mmultNxN(MMsize,ax[i][0],bx[0][j]);
    cx |= mmultNxN(MMsize,ax[i][1],bx[1][j]); // c[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j]
    if(cx!=MINIMAT0s) {
      rootc |= (1UL<<k);
      k2add_minimat(c,cx);
    }
    if(cx!=MINIMAT1s) all_ones = false;
  }
  // fix root, normalize matrix and return 
  if(rootc==NO_CHILDREN) {   // all 0s matrix is represented as an empty matrix
    assert(k2pos(c)==rootpos+1); // we wrote only root to c 
    k2setpos(c,rootpos); // delete root 
  }
  else if(all_ones && Use_all_ones_node) { // all 1s matrix is represented by the ALL_ONES root only 
    assert(rootc==ALL_CHILDREN);
    assert(k2pos(c)==rootpos+1+4*Minimat_node_ratio); // we wrote root + 4 minimats
    k2setpos(c,rootpos+1);       // discard children
    assert(k2read_node(c,rootpos)==ALL_ONES);
  }
  else {
    //root is not all_one or Use_all_ones_node is false, just write correct root node
    k2write_node(c,rootpos,rootc);
  }
}



// main entry point for matrix multiplication. 
// multiply size x size k2 compressed matrices :a and :b storing
// the result in :c, old content of :c is discarded
// :a and :b must be of size at least 2*MMsize but their content can be 
// arbitrary: all 0's, all 1's, or generic
// at exit:
//    if the result is a zero matrix: c is left empty
//    if the result is an all one's matrix: c contains a single ALL_ONES node
//    otherwise c is a node + the recursive description of its subtree  
void mmult(size_t size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c)
{
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(size>MMsize); // inputs cannot be minimats 
  assert(size%2==0);
  
  k2_free(c); // free old content and initialize as empty
  if(k2is_empty(a) ||  k2is_empty(b))
    return;  //if one matrix is all 0s the result is all 0's: nothing to be done
  
  // recursion base step
  if(size==2*MMsize) {
    mmult_base(size,a,b,c);
    return;
  }
  // case  size > 2*MMsize
  assert(size> 2*MMsize);
  node_t roota = k2read_node(a,0);
  node_t rootb = k2read_node(b,0);
  // the product of two all 1's is all 1's, but now ALL_NODES could be a pointer 
  // TODO: introduce flag for matrix using ALL_ONES
  if( (roota==ALL_ONES && a->backp==NULL) &&  (rootb==ALL_ONES && b->backp==NULL) ) {
    if(!Use_all_ones_node) //TODO: add the creation of a submatrix full of ones
      quit("Problem here: both matrices are ALL_ONES but Use_all_ones_node is false", __LINE__,__FILE__);
    assert(c->backp==NULL); // results matrices cannot have backpointers, so we can use ALL_ONES 
    k2add_node(c,ALL_ONES); // in the output matrix there are no pointers, so we can write ALL_ONES
  }  
  /*
  else if(roota==ALL_ONES && a->backp==NULL) { // further all 1s matrix optimization to be written
    left1_mmult(size,b,c);
  }
  else if(rootb==ALL_ONES && b->backp==NULL) {
    right1_mmult(size,a,c);
  }*/
  // general case: we need to split the input matrices and recurse
  else {
    split_and_rec(size,a,b,c);
  }
}


// Main entry point for right matrix vector multiplication
// multiply :asize x :asize k2 compressed matrix :a by vector :x
// storing the result in vector :y
// :a must be of size at least 2*MMsize
// :x and :y are both of :size <= :asize. 
// :a elements with index >= :size are guaranteed to be zero
// The algorithm ensures that entries of :x and :y
// with index >= :size are not accessed
// if clear_y is true y is initialized to 0, otherwise newly computed values are simply 
// added to y (this is used for example in multithread matrix-vector multiplication) 
void mvmult(size_t asize, const k2mat_t *a, size_t size, vfloat *x, vfloat *y, bool clear_y)
{
  assert(size <= asize);
  assert(asize>=2*MMsize);
  assert(asize%2==0);
  assert(a!=NULL && x!=NULL && y!=NULL);
  // initialize y to 0 if required
  if(clear_y) for(size_t i=0;i<size;i++) y[i]=0;
  if(k2is_empty(a)) return; // if a is empty the result is 0
  // call recursive multiplication algorithm based on decoding
  size_t pos = 0;
  mdecode_and_multiply(asize,a,&pos,x,y);
}


// previous recursive version of mvmult based on matrix splitting: 
// more than ten times slower than the new mvmult  
void mvmult_slow(size_t asize, const k2mat_t *a, size_t size, double *x, double *y)
{
  assert(size <= asize);
  assert(asize>=2*MMsize);
  assert(asize%2==0);
  assert(a!=NULL && x!=NULL && y!=NULL);
  // initialize y to 0
  for(size_t i=0;i<size;i++) y[i]=0;
  if(k2is_empty(a)) return; // if a is empty the result is 0
  // call recursive multiplication algorithm
  mvmult_rec(asize,a,x,y);
}



// ----------- static auxiliary functions ------------


// split input matrices and recurse matrix multiplication  
//   the input matrices are not all 0's
//   an output 0 matrix is represented as an empty matrix (no root node)   
static void split_and_rec(size_t size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c)
{
  assert(size>2*MMsize);
  assert(a!=NULL && b!=NULL && c!=NULL);
  // never called with an input empty matrix
  assert(!k2is_empty(a) && !k2is_empty(b));
  // copy *a and *b to local vars, taking care of possible back pointres
  k2mat_t atmp, btmp;
  if(a->backp!=NULL && k2read_node(a,0)==POINTER) atmp = k2jump(size,a); 
  else atmp = *a;
  if(b->backp!=NULL && k2read_node(b,0)==POINTER) btmp = k2jump(size,b); 
  else btmp = *b;
  // add subtree info if requested on-the-fly construction
  bool subinfo_added_a = false, subinfo_added_b = false;
  if(Extended_edf && a->subtinfo==NULL) {
    k2add_subtinfo_limit(size,&atmp,17);
    subinfo_added_a = true;
  }
  if(Extended_edf && b->subtinfo==NULL) {
    k2add_subtinfo_limit(size,&btmp,17);
    subinfo_added_b = true;
  }

  // split a and b into 4 blocks of half size  
  k2mat_t ax[2][2] = {{K2MAT_INITIALIZER, K2MAT_INITIALIZER}, 
                      {K2MAT_INITIALIZER, K2MAT_INITIALIZER}};
  k2mat_t bx[2][2] = {{K2MAT_INITIALIZER, K2MAT_INITIALIZER}, 
                      {K2MAT_INITIALIZER, K2MAT_INITIALIZER}}; 
  k2split_k2(size,&atmp,ax);  
  k2split_k2(size,&btmp,bx);

  // temporary matrices for products
  k2mat_t v1=K2MAT_INITIALIZER, v2=K2MAT_INITIALIZER;  
  // start building c
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_NODES as root placeholder 
  node_t rootc=NO_CHILDREN;                 // actual root node to be computed
  bool all_ones = true;
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(int k=0;k<4;k++) {
    int i=k/2, j=k%2;
    k2make_empty(&v1); k2make_empty(&v2);    // reset temporary matrices
    // a[i][0]*b[0][j]
    if(!k2is_empty(&ax[i][0]) && !k2is_empty(&bx[0][j]) ) // avoid call if one is nonzero: can be removed
      mmult(size/2, &ax[i][0], &bx[0][j], &v1);
    // a[i][1]*b[1][j]
    if(!k2is_empty(&ax[i][1]) && !k2is_empty(&bx[1][j]) )
      mmult(size/2, &ax[i][1], &bx[1][j], &v2);
    // write sum v1+v2 to c[i][j] ie block k 
    if(!k2is_empty(&v1) || !k2is_empty(&v2)) { // compute v1+v2 if at least one is nonzero
      rootc |= ((node_t )1) << k; // v1 + v2 is nonzero 
      size_t tmp = k2pos(c); // save current position to check all 1's 
      if(k2is_empty(&v2))      mcopy_full(size/2,&v1,c);  // copy v1 -> block c[i][j]  
      else if(k2is_empty(&v1)) mcopy_full(size/2,&v2,c);  // copy v2 -> block c[i][j]  
      else { // v1 and v2 are both not empty
        size_t p1=0,p2=0;
        msum_rec_plain(size/2,&v1,&p1,&v2,&p2,c);
        assert(k2pos(&v1)==p1 && k2pos(&v2)==p2); // v1 and v2 were completeley read
      }
      // check if v1+v2 is all 1's
      assert(k2pos(c)>tmp);  // since v1+v2!=0 something was written 
      if(k2read_node(c,tmp)!=ALL_ONES) all_ones = false; // v1+v2 is not all 1s
      else assert(k2pos(c)==tmp+1); // if v1+v2 is all 1s only root was written
    }
    // if v1==0 && v2==0 then v1+v2=0 and there is nothing to write and all_ones is false
    else all_ones = false;
  }
  // all computation done, deallocate temporary matrices
  k2_free(&v1);
  k2_free(&v2);
  
  // final normalization of c
  if(rootc==NO_CHILDREN) {    // case c is all 0's
    // this is an all 0 matrix, its representation is empty
    assert(k2pos(c)==rootpos+1); // only root was written to c 
    k2setpos(c,rootpos);         // delete root
  }
  else if(all_ones) {         // case c is all 1's
    assert(rootc==ALL_CHILDREN && k2pos(c) == rootpos+5);  
    // to check if this is an all 1's matrix we need to check that
    // the root is ALL_CHILDREN and that the 4 children are all 1's
    // hence are represented by a single ALL_ONES node
    for(size_t i=rootpos+1;i<rootpos+5;i++) // extra check, can be removed
       assert(k2read_node(c,i)==ALL_ONES);
    k2setpos(c,rootpos+1);      // only keep root which was already ALL_ONES
    assert(k2read_node(c,rootpos)==ALL_ONES);
  }    
  else {
    // no all 0's or all 1's, just write the correct root 
    k2write_node(c,rootpos,rootc); // fix root
  }
  // all done: deallocate subtree information if previously added
  if(subinfo_added_a) free(atmp.subtinfo);
  if(subinfo_added_b) free(btmp.subtinfo);
}

// base case of matrix vector multiplication: matrix of size 2*MMmin
// a be not all 0s (ie empty) or all 1s (these cases are handled at the previous level)
// x and y are vectors of unknown size which are never accessed
// in indices corresponding to rows/columns of a which are all 0s
// Here we call the base matrix-vector multiplication function
// using the following macro, change it to support additional sizes
#define mvmultNxN(s,a,x,y) ((s)==2 ? mvmult2x2((a),(x),(y)) : mvmult4x4((a),(x),(y)))
static void mvmult_base(size_t size, const k2mat_t *a, const vfloat *x, vfloat *y)
{
  assert(size==2*MMsize);
  assert(a!=NULL && x!=NULL && y!=NULL);
  assert(!k2is_empty(a));
  node_t roota = k2read_node(a,0);
  assert(roota!=ALL_ONES);
  // to split :a create ax[][]: the actual content will be overwritten                      
  minimat_t ax[2][2] = {{MINIMAT1s,MINIMAT1s},{MINIMAT1s,MINIMAT1s}};
  size_t posa=1; // skip the root node
  k2split_minimats(a,&posa,roota,ax);
  // split done, now matrix-vector multiply
  for(int k=0;k<4;k++) {  
    int i=k/2; int j=k%2;
    assert(size==4 || size==8); // implies that minimats are 2x2  or 4x4
    if(ax[i][j]!=MINIMAT0s) // avoid call if the block is 0
      mvmultNxN(MMsize,ax[i][j], x+j*MMsize, y + i*MMsize);
  }
}



// recursive call for matrix-vector multiplication
// compute :y += :a * :x
// :a is a non empty k2 matrix of size :size >= 2*MMsize
// :x and :y are vectors of unknown size which are never accessed
// in indices corresponding to rows/columns of :a which are all 0s
static void mvmult_rec(size_t size, const k2mat_t *a, vfloat *x, vfloat *y)
{
  assert(size>=2*MMsize);
  assert(a!=NULL && x!=NULL && y!=NULL);
  // never called with an input empty matrix
  assert(!k2is_empty(a));
  // if a is all 1s the result is easy regardless of size and the recursion stops
  node_t roota = k2read_node(a,0);
  if(roota==ALL_ONES) {
    double v = 0;
    for(size_t i=0;i<size;i++) v += x[i];
    for(size_t i=0;i<size;i++) y[i] += v;
    return;
  }
  // recursion base step
  if(size==2*MMsize)
    mvmult_base(size,a,x,y);
  else {
    // split a into 4 blocks of half size  
    k2mat_t ax[2][2] = {{K2MAT_INITIALIZER, K2MAT_INITIALIZER}, 
                        {K2MAT_INITIALIZER, K2MAT_INITIALIZER}};
    k2split_k2(size,a,ax);  
    // here we are assuming that the submatrices are in the order 00,01,10,11
    for(int k=0;k<4;k++) {
      int i=k/2, j=k%2;
      size_t z = size/2;
      if(!k2is_empty(&ax[i][j]))     // avoid call if the block is 0
        mvmult_rec(z, &ax[i][j], x+j*z, y + i*z);
    }
  }
}


// recursively decode a k2 submatrix into a list of nonzero entries
// and use each generated entry to update the matrix vector product
// Parameters:
//   size k^2 submatrix size (has the form 2^k*MMsize)
//   *c input k2mat_t structure
//   *pos position in *c where the submatrix starts
//   *x,*y input/output vector relative to the current matrix upper left corner
static void mdecode_and_multiply(size_t size, const k2mat_t *c, size_t *pos, vfloat *x, vfloat *y) {
  assert(size%2==0 && size>=2*MMsize);
  // read c root
  node_t rootc=k2read_node(c,*pos); *pos +=1;
  if(rootc==ALL_ONES) { // all 1s matrix
    double v = 0;
    for(size_t i=0;i<size;i++) v += x[i];
    for(size_t i=0;i<size;i++) y[i] += v;
    return;
  }
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(size_t k=0;k<4;k++) {  
    size_t ii = (size/2)*(k/2); size_t jj= (size/2)*(k%2);
    if(rootc & (1<<k)) { // read a submatrix
      if(size==2*MMsize) { // read a minimatrix
        minimat_t cx = k2read_minimat(c,pos);
        assert(cx!=MINIMAT0s); // should not happen
        mvmultNxN(MMsize,cx, x+jj, y+ii);
      }
      else { // decode submatrix
        mdecode_and_multiply(size/2,c,pos,x+jj, y+ii);
      }
    }
  }
}
