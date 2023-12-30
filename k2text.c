/* Routines for converting boolean matrices in text form
   to/from compressed k^2 tree representation 

   Note: internally matrix dimensions are always of the form 2^k times the size 
   of a minimatrix (those stored at the leaves of the tree), with k>0
   (somewhere this is called the k2_internal_size); the input can be of any size, 
   and the k2 matrix is padded with 0's (virtually since they are not stored)

   The conversion txt->k2 is done using an auxiliary "interleaved" array:
   each matrixc entry consits of two uint32_t (row and column indices).
   A unique entry identifier is obtained interleaving the bits of the two indices:
   as in:    r31 c31 ... r2 c2 r1 c1 r0 c0   where
   ri is the i-th bit of the row index and ci is the i-th bit of the column index
   When such interleaved values are numerically sorted the entries appear in 
   exactly the same order such entries are visited in a predorder visit of the k2 tree 
   Hence submatrices can be represented by subintervals of the interleaved array
  
   The conversion k2->txt is done using a visit of the tree in preorder
   and each time a nonzero entry is found its indices are written to the output 
   
   Currently only the minimatrix size 2 is supported
   (because of minimat_from_ia)

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "k2.h"


// prototypes of static functions
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize, size_t xsize);
static size_t mread_from_ia(uint64_t ia[], size_t n, size_t msize, k2mat_t *a);
static void mencode_ia(uint64_t *ia, size_t n, uint64_t imin, uint64_t imax, size_t size, k2mat_t *c);
static void mdecode_to_textfile(FILE *outfile, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos);
static size_t binsearch(uint64_t *ia, size_t n, uint64_t x);

// read a matrix from the text file :iname (one entry per line)
// and store it in k2 format
// the compressed matrix is stored to :a and its size to :msize
// if :xsize>0 that value is forced to be the size of k2 matrix
// return the internal size of the k2 matrix (which has the form 2**k*MMsize)
size_t mread_from_textfile(size_t *msize, k2mat_t *a, char *iname, size_t xsize)
{
  assert(iname!=NULL && a!=NULL);
  FILE *f = fopen(iname,"rt");
  if(f==NULL) quit("mread_from_file: cannot open input file",__LINE__,__FILE__);
  // generate interleaved array from input file
  size_t n; // number of entries
  uint64_t *ia = create_ia(f,&n,msize,xsize);
  assert(xsize==0 || *msize==xsize);
  fclose(f);
  // compress the matrix represented by the ia[] array into a k2mat_t structure
  size_t asize = mread_from_ia(ia,n,*msize,a);
  free(ia);
  return asize;
}



// write the content of the :size x :size k2 matrix :a to a
// text file in one entry per line format
void mwrite_to_textfile(size_t msize, size_t size, const k2mat_t *a, char *outname)
{
  assert(outname!=NULL && a!=NULL);
  assert(size>=msize);
  FILE *f = fopen(outname,"wt");
  if(f==NULL) quit("mwrite_to_file: cannot open output file",__LINE__,__FILE__);

  if(k2is_empty(a)) {  // an empty k2 matrix has no entries
    fclose(f);
    return;
  }
  size_t pos = 0;
  mdecode_to_textfile(f,msize,0,0,size,a,&pos);
  fclose(f);
  assert(pos==k2pos(a)); // check we read all the k2mat_t structure
}

// ----------- static auxiliary functions ------------

// compress the matrix of size msize represented by the interleaved
// array ia[0..n-1] into the k2mat_t structure *a 
// ia[] should be an interleaved array of length n
// the old content of :a is lost
// return the size of the k2 matrix (which has the form 2**k*MMsize)
// make sure that all entries are distinct (another option would be to 
// just remove duplicates)
static size_t mread_from_ia(uint64_t ia[], size_t n, size_t msize, k2mat_t *a)
{
  assert(ia!=NULL && a!=NULL);
  assert(n>0);                // we cannot represent an empty matrix
  assert(msize>1);
  assert( n <= msize*msize); // entries must be less that msize**2
  k2_free(a);                // free previous content of a 
  if(msize>(1<<30)) quit("mread_from_ia: matrix too large",__LINE__,__FILE__);
  size_t asize = k2get_k2size(msize);
  assert(asize>=2*MMsize);
  // count duplicates
  size_t dup=0;
  for(size_t i=1;i<n;i++)
    if(ia[i-1]==ia[i]) dup++;
  if(dup>0) {
    fprintf(stderr,"Input file contains %zu duplicate entries\n",dup);
    exit(EXIT_FAILURE);
  }
  // encode ia[] into a k2mat_t structure
  mencode_ia(ia,n,0,asize*asize,asize,a);
  return asize;
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
    assert(r==n-1);
    return n;   // replace with return r+1?
  }
  return l;
}

// recursively encode a submatrix in interleaved format  
// into a k2mat_t structure
// Parameters:
//   ia[0,n-1] array containing the distinct interleaved entries 
//   imin, imax range of values in ia[0,n-1]
//   size  submatrix size (has the form 2^k*MMsize)
//   *c output k2mat_t structure to be filled in dfs order 
static void mencode_ia(uint64_t *ia, size_t n, uint64_t imin, uint64_t imax, size_t size, k2mat_t *c) {
  assert(ia!=NULL);
  assert(n>0);
  assert(imin>=0 && imax>=0);
  assert(imin<imax);
  assert(ia[0]>=imin); 
  assert(ia[n-1]<imax);
  assert(size%2==0 && size>=2*MMsize);
  //\\ printf("Size=%zu, n=%zu, imin=%ld imax=%ld\n",size,n,imin,imax);
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
  int64_t mid = left+range;
  int64_t right = imax - range;
  assert(right==mid+range);
  // determine range in ia[] of the 4 submatrices entries
  size_t imid = binsearch(ia,n,mid);    //   first entry of A[10]
  size_t ileft = imid>0 ? binsearch(ia,imid,left):0; // first entry of A[01]
  size_t iright = imid<n? binsearch(ia+imid,n-imid,right)+imid:n; // first entry of A[11]
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



// interleaves two 32 bits integers in a single uint64_t 
// the bits of a (row index) are more significant than
// those of b (column index) because of how we number submatrices
static uint64_t bits_interleave(int64_t a, int64_t b)
{
  uint64_t r = 0;
  assert(a<=UINT32_MAX && b <= UINT32_MAX);
  int c = 0;
  while(a!=0 || b!=0) {
    r |= (b&1)<<c++;
    r |= (a&1)<<c++;
    a >>= 1; b>>=1;  
    assert(c<=64);
  }
  return r;
}

static int uint64_cmp(const void *p, const void *q)
{
  const uint64_t *a = p;
  const uint64_t *b = q;
  
  if(*a < *b) return -1;
  else if(*a > *b) return 1;
  return 0;
}

// create and return interleaved array from the list of entries in a text file
// the matrix size stored in :msize is computed as follows: 
//  if xsize==0 *msize = largest index + 1
//  if xsize>0 that value is forced to be the matrix size (all indexes must be <xsize)
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize, size_t xsize)
{
  uint32_t maxentry = 0; // largest entry in the file
  size_t size=10;      // current size of ia[]
  size_t i=0;          // elements in ia[]
  uint64_t *ia = malloc(size*sizeof(*ia));
  if(ia==NULL) quit("create_ia: malloc failed",__LINE__,__FILE__);
    
  int64_t a,b; size_t line=0;  
  while(true) {
    line++;
    int e = fscanf(f,"%" SCNd64 " %" SCNd64,&a,&b);
    if(e==EOF) break;
    // check input
    if(e!=2) {
      fprintf(stderr,"Invalid file content at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(a<0 || b<0) {
      fprintf(stderr,"Negative index at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(a>UINT32_MAX || b>UINT32_MAX) {
      fprintf(stderr,"Index too large at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(xsize>0 && (a>=xsize || b>=xsize)) {
      fprintf(stderr,"Index larger than the assigned size at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    // update maxentry
    if(a>maxentry) maxentry=a;
    if(b>maxentry) maxentry=b;
    // compute interleaved value
    uint64_t entry = bits_interleave(a,b);
    // enlarge ia if necessary
    if(i==size) {
        size = size*2;
        ia = realloc(ia,size*sizeof(*ia));
        if(ia==NULL) quit("create_ia: realloc failed",__LINE__,__FILE__);
    }
    assert(size>i);
    ia[i++] = entry;
  }
  // final resize
  size = i;
  ia = realloc(ia,size*sizeof(*ia));
  if(ia==NULL) quit("create_ia: realloc failed",__LINE__,__FILE__);
  // sort interleaved entries
  qsort(ia, size, sizeof(*ia), &uint64_cmp);
  //\\ for(i=0;i<size;i++) printf("%ld ", ia[i]); puts("");
  // save output parameters   
  if(xsize==0) { // if xsize==0 size is largest index + 1
    if(maxentry+1>SIZE_MAX)  // highly unlikely, but you never know... 
      quit("create_ia: cannot represent matrix size",__LINE__,__FILE__);
    *msize = maxentry+1;
  }
  else {  // if parameter xsize>0 that is the desired matrix size
    assert(maxentry<xsize);
    *msize = xsize;
  }
  *n = size;
  return ia;  
}


// -----------------------------

// recursively decode a k2 submatrix into a list of entries
// written to a text file
// Parameters:
//   f output file 
//   msize actual file of the matrix
//   i,j submatrix top left corner
//   size k^2 submatrix size (has the form 2^k*MMsize)
//   *c input k2mat_t structure
//   *pos position in *c where the submatrix starts
static void mdecode_to_textfile(FILE *outfile, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos) {
  assert(size%2==0 && size>=2*MMsize);
  assert(i%MMsize==0 && j%MMsize==0);
  assert(i>=0 && j>=0 && i<msize+2*size && j<msize+2*size);
  // read c root
  node_t rootc=k2read_node(c,*pos); *pos +=1;
  if(rootc==ALL_ONES) { // all 1s matrix
    for(size_t ii=0; ii<size; ii++)
      for(size_t jj=0; jj<size; jj++)
        if(i+ii<msize && j+jj<msize) { 
          int e = fprintf(outfile,"%zu %zu\n",i+ii,j+jj);
          if(e<0) quit("mdecode_to_textfile: fprintf failed",__LINE__,__FILE__);
        }
    return;
  }
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(int k=0;k<4;k++) {  
    size_t ii = i + (size/2)*(k/2); size_t jj= j + (size/2)*(k%2);
    if(rootc & (1<<k)) { // read a submatrix
      if(size==2*MMsize) { // read a minimatrix
        minimat_t cx = k2read_minimat(c,pos);
        assert(cx!=MINIMAT0s); // should not happen
        minimat_to_text(outfile,msize,ii,jj,size/2,cx);
      }
      else { // decode submatrix
        mdecode_to_textfile(outfile,msize,ii,jj,size/2,c,pos);
      }
    }
  }
}

