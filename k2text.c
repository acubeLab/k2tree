/* Routines for converting binary matrices in text form, ie 
   a list of directed arcs represented as a pair of vertices,
   to the compressed k^2 tree representation 

   Matrix dimensions are assumed to be power of 2, of size at least
   2*MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "k2aux.c"
#include "bbm.h"

// prototypes of static functions
static void mencode_ia(uint64_t *ia, size_t n, int64_t imin, int64_t imax, int size, k2mat_t *c);
static size_t binsearch(uint64_t *ia, size_t n, uint64_t x);


// write the content of the :size x :size k2 matrix :a to the bbm matrix :m 
// of size msize*msize. It is assumed m was already correctly initialized and allocated
void mwrite_to_bbm(uint8_t *m, int msize, int size, const k2mat_t *a)
{
  assert(size>=msize);
  if(k2is_empty(a)) {  // an empty k2 matrix is all 0s
    byte_to_bbm(m,msize,0,0,msize,0); // fill m with 0s
    return;
  }
  byte_to_bbm(m,msize,0,0,msize,2); // fill m with illegal value 2
  size_t pos = 0;
  mdecode(m,msize,0,0,size,a,&pos);
  assert(pos==k2pos(a)); // check we read all the k2mat_t structure
}


// compress the matrix *m of size msize represented by the interleaved
// array ia[0..n-1] into the k2mat_t structure *a 
// ia[] should be an interleaved array of length n
// the old content of :a is lost
// return the size of the k2 matrix (which has the form 2**k*MMsize)
int mread_from_ia(uint64_t ia[], size_t n, int msize, k2mat_t *a)
{
  assert(ia!=NULL && a!=NULL);
  assert(n>0);          // we cannot represent an empty matrix
  assert(msize>1);
  assert( n <= ((size_t)msize)*msize); // we cannot represent a matrix larger than MMsize
  k2_free(a); // free previous content of a 
  if(msize>(1<<30)) quit("mread_from_ia: matrix too large",__LINE__,__FILE__);
  int asize = k2get_k2size(msize);
  assert(asize>=2*MMsize);
  int64_t asize2 = ((int64_t) asize) *asize;
  // encode ia[] into a k2mat_t structure
  mencode_ia(ia,n,0,asize2,asize,a);
  return asize;
}


// ----------- static auxiliary functions ------------

// recursively encode a submatrix in interleaved format  
// into a k2mat_t structure
// Parameters:
//   ia[0,n-1] array containing the interleaved arcs 
//   imin, imax range of values in ia[0,n-1]
//   size  submatrix size (has the form 2^k*MMsize)
//   *c output k2mat_t structure to be filled in dfs order 
static void mencode_ia(uint64_t *ia, size_t n, int64_t imin, int64_t imax, int size, k2mat_t *c) {
  assert(n>0);
  assert(imin>=0 && imax>=0);
  assert(imin<imax);
  assert(ia[0]>=imin && ia[n-1]<imax);
  assert(size%2==0 && size>=2*MMsize);
  // case of a full submatrix
  if(n==size*size) {
    k2add_node(c,ALL_ONES);   // submatrix is full 
    return;
  }
  // determine range of submatrices 
  int64_t range = imax-imin;
  assert(range%4==0);
  range = range/4;
  int64_t left = imin + range;
  int64_t mid = imin
  int64_t right = imax - range;
  // determine range in ia[] of the 4 submatrices entries
  size_t imid = binsearch(ia,n,mid);    //   first entry of A[10]
  size_t ileft = binsearch(ia,imid,left); // first entry of A[01]
  size_t iright = binsearch(ia+imid,n-imid,right); // first entry of A[11]
  // the four submatrices are: 
  //    ia[0,ileft-1], ia[ileft,imid-1], ia[imid,iright-1], ia[iright,n-1]
  // and contain values in the ranges
  //    [imin,left), [left,mid), [mid,right), [right,imax)
  // start building c
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_ONES as root placeholder 
  node_t rootc=NO_CHILDREN;                 // actual root node to be computed
  // here we are assuming that the submatrices are in the order 00,01,10,11

  if(ileft>0) { // submatrix 00 is not empty
    rootc |= (1<<0); // set 00 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia,ileft,imin,left,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia,ileft,imin,left,size/2);
      k2add_minimat(c,cx);
    }
  }

  if(ileft<imid) { // submatrix 01 is not empty
    rootc |= (1<<1); // set 01 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia+ileft,imid-ileft,left,mid,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia+ileft,imid-ileft,left,mid,size/2);
      k2add_minimat(c,cx);
    }
  }

  if(iright>imid) { // submatrix 10 is not empty
    rootc |= (1<<2); // set 10 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia+imid,iright-imid,mid,right,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia+imid,iright-imid,mid,right,size/2);
      k2add_minimat(c,cx);
    }
  }

  if(iright<n) { // submatrix 11 is not empty
    rootc |= (1<<3); // set 11 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia+iright,n-iright,right,imax,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia+iright,n-iright,right,imax,size/2);
      k2add_minimat(c,cx);
    }
  }
  assert(rootc!=NO_CHILDREN); // at least one submatrix is not empty
  k2write_node(c,rootpos,rootc); // fix root 
}


// in a sorted uint64_t array ia[0,n-1] containing distinct values find 
// the first entry >= x using binary search 
// return n if no such entry exists
static size_t binsearch(uint64_t *ia, size_t n, uint64_t x) {
  assert(ia!=NULL && n>0);
  size_t l=0, r=n-1;
  while(l<r) {
    size_t m = (l+r)/2;
    if(ia[m]<x) l=m+1;
    else if(ia[m]==x) return m;
    else r=m; // replace with r = m-1 and later return r+1?
  }
  assert(l==r);
  if(ia[l]<x) {
    assert(r==n-1)
    return n;   // replace with return r+1?
  }
  return l;
}


// recursively decode a k2 submatrix into a binary submatrix
// m[i,i+size)[j,j+size) in one-byte per entry format into 
// Parameters:
//   m output matrix of (overall) size mszie*msize
//   i,j submatrix top left corner
//   size submatrix size (has the form 2^k*MMsize)
//   *c input k2mat_t structure
//   *pos position in *c where the submatrix starts
static void mdecode(uint8_t *m, int msize, int i, int j, int size, const k2mat_t *c, size_t *pos) {
  assert(size%2==0 && size>=2*MMsize);
  assert(i%MMsize==0 && j%MMsize==0);
  assert(i>=0 && j>=0 && i<msize+2*size && j<msize+2*size);
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
      if(k2is_empty(&v2))      mcopy(size/2,&v1,c);  // copy v1 -> block c[i][j]  
      else if(k2is_empty(&v1)) mcopy(size/2,&v2,c);  // copy v2 -> block c[i][j]  
      else { // v1 and v2 are both not empty
        size_t p1=0,p2=0;
        msum_rec(size/2,&v1,&p1,&v2,&p2,c);
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
}
