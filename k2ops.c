/* Routines for arithmetic/logical operations on binary matrices 
   represented as k^2 trees 

   This file contains the defintion of the operations that make use of
   the basic operations defined in k2aux.c

   Matrix dimensions are assumed to be power of 2, of size at least
   2*MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "k2aux.c"
#include "bbm.h"

static void mdecode(uint8_t *m, int msize, int i, int j, int size, k2mat_t *c, size_t *pos);
static void mencode(uint8_t *m, int msize, int i, int j, int size, k2mat_t *c);


// write the content of the :size x :size k2 matrix :a to the bbm matrix :m 
// of size msize*msize. It is assumed m was already correctly initialized 
void mwrite_to_bbm(uint8_t *m, int msize, int size, k2mat_t *a)
{
  assert(size>=msize);
  byte_to_bbm(m,msize,0,0,msize,0); // fill m with illegal value 2
  size_t pos = 0;
  mdecode(m,msize,0,0,size,a,&pos);
  assert(pos==k2pos(a)); // check we read all the k2mat_t structure
}


// read the uncompressed matrix *m of size msize into the k2mat_t structure *a 
// m should be an integer array of size msize*msize 
// return the size of the k2 matrix (which has the form 2**k*MMsize)
int mread_from_bbm(uint8_t *m, int msize, k2mat_t *a)
{
  assert(msize>1);
  assert(k2is_empty(a));
  assert(!a->read_only);
  if(msize>(1<<30)) quit("k2read_uncompressed: matrix too large",__LINE__,__FILE__);
  int size = k2get_k2size(msize);
  assert(size>=2*MMsize);
  // read matrix m into a
  mencode(m,msize,0,0,size,a);
  return size;
}

// copy matrix a: to b: used instead of sum when one of the matrices is all 0s
void mcopy(int size, const k2mat_t *a, k2mat_t *b)
{
  assert(size>MMsize);  // required by k2copy_rec 
  assert(a!=NULL && b!=NULL);
  assert(!k2is_empty(a));
  assert(!b->read_only);
  size_t pos = 0;
  k2copy_rec(size,a,&pos,b);
 }


// NO LONGER USED: base case solved in the general case
// base case of matrix sum: matrices of size 2*MMmin
// a and b must not be 0s (ie empty)
//   (if a or b is all 0s this function is not called because the sum is a copy) 
// a and b can be all 1's 
// the output matrix c is normalized as ususal
//  if c is all 0s nothing is written (this should not happen because the sum is an OR)
//  if c is all 1s just the root ALL_ONES is written
//  otherwise we follow the standard rule: root node + nonzero minisize matrices   
void _xxx_msum_base(int size, const k2mat_t *a, size_t *posa, 
                         const k2mat_t *b, size_t *posb, k2mat_t *c)
{
  assert(size==2*MMsize);
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(!k2is_empty(a) && !k2is_empty(b));
  // take care of all 1s matrices: read root without advancing positions
  node_t roota = k2read_node(a,*posa);
  node_t rootb = k2read_node(b,*posb);
  size_t tmp1,tmp2; // not used, just to pass to k2dfs_visit
  if(roota==ALL_ONES) {
    k2add_node(c,ALL_ONES); 
    *posa+=1; // reach the end of a
    k2dfs_visit(size,b,posb,&tmp1,&tmp2); //reach end of b but ignore its content
    return;
  }
  else if(rootb==ALL_ONES) {
    k2add_node(c,ALL_ONES); 
    *posb+=1; // same as above with a and b swapped
    k2dfs_visit(size,a,posa,&tmp1,&tmp2); //scan but ignore a content
    return;
  }
  // a and b are not all 1s, merge children
  *posa += 1; *posb += 1; // skip root nodes already read
  node_t rootc= roota | rootb; // root node of c, correct except when c is all 1s
  size_t rootpos = k2add_node(c,rootc);  // save root and its position
  minimat_t ax, bx;
  bool all_ones=true; // true if all submatrices cx[i][j] are all 1's
  for(int k=0;k<4;k++) {
    if(roota | (1<<k)) {
      ax = k2read_minimat(a,*posa); *posa += Minimat_node_ratio;
    } else ax = MINIMAT0s;
    if(rootb | (1<<k)) {
      bx = k2read_minimat(b,*posb); *posb += Minimat_node_ratio;
    } else bx = MINIMAT0s;
    minimat_t cx = ax | bx;   // compute or of correposnding children
    if(cx!=MINIMAT0s) {       // save it if nonzero
      rootc |= (1<<k);
      k2add_minimat(c,cx);
    }
    if(cx!=MINIMAT1s) all_ones = false;
  }
  // possible normalization if c is all 1s
  assert(rootc!=NO_CHILDREN); // at least one child is nonzero
  if(all_ones) {
    assert(rootc==ALL_CHILDREN);
    assert(k2pos(c)==rootpos+1+4*Minimat_node_ratio); // we wrote root + 4 minimat
    k2setpos(c,rootpos);       // discard current root and children
    k2write_node(c,rootpos,ALL_ONES); // write ALL_ONES as root
  }
  return; 
}


// recursive sum of two k2 (sub)matrices 
// a and b must not be all 0s (sub)matrices
//   (if a or b is all 0s this function is not called because the sum is a copy) 
// a and b can be all 1's 
// the output matrix c is normalized as ususal
//  if c is all 0s nothing is written (this should not happen because the sum is an OR)
//  if c is all 1s just the root ALL_ONES is written
//  otherwise we follow the standard rule: root node + subtrees in DFS order   
void msum_rec(int size, const k2mat_t *a, size_t *posa, 
                         const k2mat_t *b, size_t *posb, k2mat_t *c)
{
  assert(size>=2*MMsize); // there must be the root node 
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(!k2submatrix_empty(a,*posa) && !k2submatrix_empty(b,*posb));
  // take care of all 1s matrices: read root without advancing positions
  node_t roota = k2read_node(a,*posa);
  node_t rootb = k2read_node(b,*posb);
  size_t tmp1=0,tmp2=0; // not used, defined to pass to k2dfs_visit
  if(roota==ALL_ONES) {
    k2add_node(c,ALL_ONES); 
    *posa+=1; // reach the end of a
    k2dfs_visit(size,b,posb,&tmp1,&tmp2); //reach end of b but ignore its content
    return;
  }
  else if(rootb==ALL_ONES) {
    k2add_node(c,ALL_ONES); 
    *posb+=1; // same as above with a and b swapped
    k2dfs_visit(size,a,posa,&tmp1,&tmp2); //scan but ignore a content
    return;
  }
  assert(roota!=ALL_ONES && rootb!=ALL_ONES);
  // a and b are not all 1s, merge children
  *posa += 1; *posb += 1; // skip root nodes already read
  node_t rootc = roota | rootb; // root node of c, correct except when c is all 1s
  assert(rootc!=NO_CHILDREN);   // at least one child is nonzero
  size_t rootpos = k2add_node(c,rootc);  // save c root and its position
  bool all_ones=true;  // true if all submatrices cx[i][j] are all 1's
  if(size==2*MMsize) { // children are minimat matrices
    minimat_t ax[2][2], bx[2][2];
    k2split_minimats(a,posa,roota,ax);
    k2split_minimats(b,posb,rootb,bx);
    for(int k=0;k<4;k++) {
      minimat_t cx = ax[k/2][k%2] | bx[k/2][k%2]; // compute bitwise or of corresponding minimat
      if (cx != MINIMAT0s)
      { // save cx if nonzero
        assert(rootc & (1 << k)); // implied by  cx = ax | bx 
        k2add_minimat(c, cx);
      }
      else assert((rootc & (1 << k)) == 0);
      if (cx != MINIMAT1s) all_ones = false;
    }
    assert(all_ones==false || rootc==ALL_CHILDREN); // all_ones => rootc==ALL_CHILDREN
    if(all_ones) assert(k2pos(c)==rootpos+1+4*Minimat_node_ratio);
  }
  else { // size>2*MMsize: children are k2 matrices, possibly use recursion
    // we could split a and b in 4 submatrices and sum them, but it is more efficient 
    // to compute the sum without building submatrices (which requires a scan of a and b)
    for(int k=0;k<4;k++) {
      if (roota | (1 << k)) {
        if(rootb | (1 << k)) 
          msum_rec(size/2,a,posa,b,posb,c); // k-th child of c is sum of kth children of a and b
        else 
          k2copy_rec(size/2,a,posa,c); // k-th child of c is kth child of a
      }
      else if (rootb | (1 << k))
        k2copy_rec(size/2,b,posb,c); // k-th child of c is kth child of b
      else 
        assert( (rootc & (1 << k)) == 0); // both children are 0s, nothing to do 
    }
    // check if all children of rootc are all 1's
    // done checking that the subtree contains just 5 nodes 
    all_ones =  rootc==ALL_CHILDREN && k2pos(c)==rootpos+5;
  }
  // normalize if c is all 1s (regardless of size)
  if(all_ones) {
    assert(rootc==ALL_CHILDREN);
    k2setpos(c,rootpos);       // discard current root and children
    k2write_node(c,rootpos,ALL_ONES); // write ALL_ONES as root
  }
  return; 
}

// main entry point for matrix addition
// sum size x size k2 compressed matrices :a and :b storing
// the result in c, must be called with an initialized and empty :c
// :a and :b must be of size at least 2*MMsize but their content can be 
// arbitrary: all 0's, all 1's, or generic
// at exit:
//    if the result is a zero matrix c is left empty
//    if the result is an all one's matrix c contains a single ALL_ONES node
//    otherwise c is a node + the recursive description of its subtree  
void msum(int size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c)
{
  assert(size>=2*MMsize);
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(k2is_empty(c) && !c->read_only);

  if(k2is_empty(a) && k2is_empty(b)) 
    return;                 // if a==0 && b==0: c=0
  else if(k2is_empty(b))      
    mcopy(size,a,c);        // if b==0: c=a
  else if(k2is_empty(a))    
    mcopy(size,b,c);        // if a==0: c=b
  else { // a and b are both not empty, call msum_rec
    size_t posa=0,posb=0;
    msum_rec(size,a,&posa,b,&posb,c);
    assert(posa==k2pos(a) && posb==k2pos(b)); // a and b were completeley read
    assert(!k2is_empty(c));  // implied by a+b!=0
  }
  return;
}

// recursive test for equality test of two k2 matrices both nonzero
// see mequals for details
int mequals_rec(int size, const k2mat_t *a, size_t *posa, 
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
// these two values differ by one (these can be seen in the returned value)
// :a and :b must be of size at least 2*MMsize but their content can be
// arbitrary: all 0's, all 1's, or generic
// note: here all 0's matrices are considered of depth 1 even if they are empty
int mequals(int size, const k2mat_t *a, const k2mat_t *b)
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


// base case of matrix multiplication: matrices of size 2*MMmin
// a and b must be not all 0s (ie empty)
//   (if a or b is 0s this function is not called because the product is 0) 
// a and b can be all 1's 
// the output matrix c is normalized as usual:
//  if c is all 0s nothing is written
//  if c is all 1s the root ALL_ONES is written
//  otherwise we follow the standard rule: root node + nonzero minisize matrices  
void mmult_base(int size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c)
{
  assert(size==2*MMsize);
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(!k2is_empty(a) && !k2is_empty(b));
  assert(!c->read_only);
  // c is always an empty matrix because a product is never written directly to a result matrix 
  assert(k2pos(c)==0);  
  // initialize a[][] and b[][] to cover the case when the matrix is all 1s                       
  minimat_t ax[2][2] = {{MINIMAT1s,MINIMAT1s},{MINIMAT1s,MINIMAT1s}};
  minimat_t bx[2][2] = {{MINIMAT1s,MINIMAT1s},{MINIMAT1s,MINIMAT1s}};
  node_t roota = k2read_node(a,0);
  node_t rootb = k2read_node(b,0);
  //both matrices are all 1s ?
  if(roota==ALL_ONES && rootb == ALL_ONES) {
    k2add_node(c,ALL_ONES);
    return;
  }
  // ??? here possible code for case when one matrix is all 1s and the other is not
  // split a and b
  size_t posa=1,posb=1; // we have already read the root node
  if(roota!=ALL_ONES)   // case ALL_ONES is coverede by initialization above
    k2split_minimats(a,&posa,roota,ax);
  if(rootb!=ALL_ONES) 
    k2split_minimats(b,&posb,rootb,bx);
  // split done, now multiply and store c 
  // optimization on all 1's submatrices still missing 
  minimat_t cx[2][2];
  bool all_ones=true; // true if all submatrices cx[i][j] are all 1's
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_NODES as root placeholder 
  node_t rootc=NO_CHILDREN;        // actual root node to be computed
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(int k=0;k<4;k++) {  
    int i=k/2; int j=k%2;
    assert(size==4); // implies that minimats are 2x2 
    cx[i][j]  = mmult2x2(ax[i][0],bx[0][j]);
    cx[i][j] |= mmult2x2(ax[i][1],bx[1][j]);
    if(cx[i][j]!=MINIMAT0s) {
      rootc |= (1<<k);
      k2add_minimat(c,cx[i][j]);
    }
    if(cx[i][j]!=MINIMAT1s) all_ones = false;
  }
  // fix root normalize matrix and return 
  if(rootc==NO_CHILDREN) {   // all 0s matrix is represented as a an empty matrix
    assert(k2pos(c)==rootpos+1); // we wrote only root to c 
    k2setpos(c,rootpos); // delete root 
  }
  else if(all_ones) {    // all 1s matrix is represented by the ALL_ONES root only 
    assert(rootc==ALL_CHILDREN);
    assert(k2pos(c)==rootpos+1+4*Minimat_node_ratio); // we wrote root + 4 minimats
    k2setpos(c,rootpos+1);       // discard children
    assert(k2read_node(c,rootpos)==ALL_ONES);
  }
  else k2write_node(c,rootpos,rootc); // just fix root 
}


// split input matrices and recurse matrix multiplication  
//   an input 0 matrix is represented as NULL
//   an output 0 matrix is represented as an empty matrix (no root node)   
static void split_and_rec(int size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c)
{
  assert(size>2*MMsize);
  assert(a!=NULL && b!=NULL && c!=NULL);
  // never called with an input empty matrix
  assert(!k2is_empty(a) && !k2is_empty(b));
  // split a and b into 4 blocks of half size  
  k2mat_t ax[2][2] = {{K2MAT_INITIALIZER, K2MAT_INITIALIZER}, 
                      {K2MAT_INITIALIZER, K2MAT_INITIALIZER}};
  k2mat_t bx[2][2] = {{K2MAT_INITIALIZER, K2MAT_INITIALIZER}, 
                      {K2MAT_INITIALIZER, K2MAT_INITIALIZER}}; 
  k2split_k2(size,a,ax);  
  k2split_k2(size,b,bx);
  // temporary matrices for products
  k2mat_t v1=K2MAT_INITIALIZER, v2=K2MAT_INITIALIZER;
  
  // start building c
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_NODES as root placeholder 
  node_t rootc=NO_CHILDREN;                 // actual root node to be computed
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
    // sum to c[i][j] ie block k 
    if(!k2is_empty(&v1) || !k2is_empty(&v2)) { // sum if at least one is nonzero
      rootc |= ((node_t )1) << k; // v1 + v2 is nonzero 
      // compute v1+v2
      if(k2is_empty(&v2)) 
        mcopy(size,&v1,c);  // copy v1 -> block c[i][j]  
      else if(k2is_empty(&v1)) 
        mcopy(size,&v2,c);  // copy v2 -> block c[i][j]  
      else { // v1 and v2 are both not empty
        size_t tmp = k2pos(c),p1=0,p2=0;
        msum_rec(size,&v1,&p1,&v2,&p2,c);
        assert(k2pos(&v1)==p1 && k2pos(&v2)==p2); // v1 and v2 were completeley read
        assert(tmp<k2pos(c));  // implied by v1+v2!=0
      }
    }
    // if v1==0 && v2==0 then v1+v2=0 and there is nothing to write 
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
  // case c is all 1's
  else if(rootc==ALL_CHILDREN && k2pos(c) == rootpos+5)  {
    // to check if this is an all 1's matrix we need to check that
    // the root is ALL_CHILDREN and that the 4 children are all 1's
    // hence are represented by a single ALL_ONES node
    for(size_t i=rootpos+1;i<rootpos+5;i++)
       assert(k2read_node(c,i)==ALL_ONES);
    k2setpos(c,rootpos+1);      // only keep root which was already ALL_ONES
    assert(k2read_node(c,rootpos)==ALL_ONES);
  }    
  else {
    // no all 0's or all 1's, just write the correct root 
    k2write_node(c,rootpos,rootc); // fix root
  }  
}



// main entry point for matrix multiplication. 
// multiply size x size k2 compressed matrices :a and :b storing
// the result in c, must be called with an initialized and empty :c
// :a and :b must be of size at least 2*MMsize but their content can be 
// arbitrary: all 0's, all 1's, or generic
// at exit:
//    if the result is a zero matrix c is left empty
//    if the result is an all one's matrix c contains a single ALL_ONES node
//    otherwise c is a node + the recursive description of its subtree  
void mmult(int size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c)
{
  assert(a!=NULL && b!=NULL && c!=NULL);
  assert(size>MMsize); // inputs cannot be minimats 
  assert(size%2==0);
  assert(k2is_empty(c));
  
  //if one matrix is all 0s (empty) the result is all 0's: nothing to be done  
  if(k2is_empty(a) ||  k2is_empty(b))
    return;
  
  // recursion base step
  if(size==2*MMsize) {
    mmult_base(size,a,b,c);
    return;
  }
  // general case 
  assert(size> 2*MMsize);
  node_t roota = k2read_node(a,0);
  node_t rootb = k2read_node(b,0);
  // the product of two all 1's is all 1's 
  if(roota==ALL_ONES &&  rootb==ALL_ONES) {
    k2add_node(c,ALL_ONES);
  }
  // further all 1s matrix optimization to be written
  /*
  else if(na==ALL_ONES) {
    left1_mmult(size,b,c);
  }
  else if(nb==ALL_ONES) {
    right1_mmult(size,a,c);
  }*/
  split_and_rec(size,a,b,c);
}

// return statistics on matrix a
// write number of used pos,nomde, and minimats in the variables passed by reference
// and return the number of levels
int mstats(int size, const k2mat_t *a, size_t *pos, size_t *nodes, size_t *minimats)
{
  *pos=*nodes=*minimats=0;
  if(!k2is_empty(a)) k2dfs_visit(size,a,pos,nodes,minimats);
  int eq = mequals(size,a,a);
  assert(eq<0);
  return -eq;
}


// recursively encode a binary submatrix m[i,i+size)[j,j+size) given in one-byte 
// per entry format into a k2mat_t structure
// Parameters:
//   m matrix to encode of size mszie*msize
//   i,j submatrix top left corner
//   size submatrix size (has the form 2^k*MMsize)
//   *c output k2mat_t structure to be filled in dfs order 
static void mencode(uint8_t *m, int msize, int i, int j, int size, k2mat_t *c) {
  assert(size%2==0 && size>=2*MMsize);
  assert(i%MMsize==0 && j%MMsize==0);
  assert(i>=0 && j>=0 && i<msize+2*size && j<msize+2*size);
  // if we are outside m it's an all 0 submatrix and there is nothing to do
  if(i>=msize || j>=msize) return;
  // start building c
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_NODES as root placeholder 
  node_t rootc=NO_CHILDREN;                 // actual root node to be computed
  bool all_ones=true;  // true if all submatrices cx[i][j] are all 1's
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(int k=0;k<4;k++) {  
    int ii = i + (size/2)*(k/2); int jj= j + (size/2)*(k%2);
    if(size==2*MMsize) {
      minimat_t cx = minimat_from_bbm(m,msize,ii,jj,size/2);
      if(cx!=MINIMAT0s) {
        rootc |= (1<<k);
        k2add_minimat(c,cx);
      }
      if(cx!=MINIMAT1s) all_ones = false;
    } // size>2*MMsize: use recursion
    else {
      size_t tmp = k2pos(c); // save current position
      mencode(m,msize,ii,jj,size/2,c);
      if(tmp!=k2pos(c)) { // something was written
        assert(k2pos(c)>tmp);
        rootc |= (1<<k);
        if(k2read_node(c,tmp)!=ALL_ONES) all_ones = false;
        else assert(k2pos(c)==tmp+1);
      }
      else all_ones = false; // nothing was written
    }
  }
  // compute (or check) all_ones
  if(size!=2*MMsize) {
    if(all_ones)  assert(rootc==ALL_CHILDREN && k2pos(c)==rootpos+5);}
  else {
    if(all_ones) assert(k2pos(c)==rootpos+1+4*Minimat_node_ratio); }
  // fix root and normalize matrix
  if(rootc==NO_CHILDREN) { // all 0s matrix is represented as an empty matrix
    assert(k2pos(c)==rootpos+1); // we wrote only root to c 
    k2setpos(c,rootpos); // delete root  
  }
  else if(all_ones) {   // all 1s matrix is represented by the ALL_ONES root only
    assert(rootc==ALL_CHILDREN);
    k2setpos(c,rootpos+1);       // discard children
    assert(k2read_node(c,rootpos)==ALL_ONES); // ALL_ONES was the defualt root 
  }
  else k2write_node(c,rootpos,rootc); // just fix root 
}

// recursively decode a k2 submatrix into a binary submatrix
// m[i,i+size)[j,j+size) in one-byte per entry format into 
// Parameters:
//   m output matrix of (overall) size mszie*msize
//   i,j submatrix top left corner
//   size submatrix size (has the form 2^k*MMsize)
//   *c input k2mat_t structure
//   *pos position in *c where the submatrix starts
static void mdecode(uint8_t *m, int msize, int i, int j, int size, k2mat_t *c, size_t *pos) {
  assert(size%2==0 && size>=2*MMsize);
  assert(i%MMsize==0 && j%MMsize==0);
  assert(i>=0 && j>=0 && i<msize+2*size && j<msize+2*size);
  printf("mdecode: i=%d j=%d size=%d\n",i,j,size);
  // read c root
  node_t rootc=k2read_node(c,*pos); *pos +=1;
  if(rootc==ALL_ONES) { // all 1s matrix
    byte_to_bbm(m,msize,i,j,size,1);
    return;
  }
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(int k=0;k<4;k++) {  
    int ii = i + (size/2)*(k/2); int jj= j + (size/2)*(k%2);
    if(rootc & (1<<k)) { // read a submatrix
      if(size==2*MMsize) { // read a minimatrix
        minimat_t cx = k2read_minimat(c,*pos); *pos += Minimat_node_ratio;
        assert(cx!=MINIMAT0s); // should not happen
        minimat_to_bbm(m,msize,ii,jj,size/2,cx);
      }
      else { // decode submatrix
        mdecode(m,msize,ii,jj,size/2,c,pos);
      }
    }
    else { // the k-th chilldren is 0: write a 0 submatrix
      // if m initialized with 0s we can skip the following line and save some writes
      byte_to_bbm(m,msize,ii,jj,size/2,0);
    }
  }
}
