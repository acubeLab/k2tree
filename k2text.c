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
   
   Note that the size of the interleaved array is equal to the number of 
   nonzeros, which we assume is less than 2^64, hence indices in the
   array can be stored in a size_t. However, the single entry store the 
   row and column index so it must be able to store a number of bits equal 
   to (2 x bits in a single index).
   Currently the maximum allowed size is 2^32, so each index takes 32 bits
   and the interleaved array can be of int64_t's. To support larger 
   matrices, say up to 2^40, the entries of the interleaved array ia[]
   and the related variables (imin,left,mid,right) must be enlarged.
   This can be done using uint128_t for the scalars and an appropriate 
   byte array for ia[].
   
   Recall than when working with values >= 2^32 stored in an uint64_t 
   we cannot safely compute products: this is why we have the functions
   altb2 and aeqb2 testing whether a<b*b or a==b*b without multiplications 
    
   The conversion k2->txt is done using a visit of the tree in preorder
   and each time a nonzero entry is found its indices are written to the output 
   
   Currently only the minimatrix sizes 2 and 4 are supported

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <inttypes.h>
#include "k2.h"


// prototypes of static functions
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize, size_t xsize);
static size_t mread_from_ia(uint64_t ia[], size_t n, size_t msize, k2mat_t *a);
static void mencode_ia(uint64_t *ia, size_t n, uint64_t imin, size_t size, k2mat_t *c);
static void mdecode_to_textfile(FILE *outfile, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos);
static size_t binsearch(uint64_t *ia, size_t n, uint64_t x);
static __inline__ bool a_eq_b2(uint64_t a, uint64_t b);
static __inline__ bool a_lt_b2(uint64_t a, uint64_t b);



// read a matrix from the text file :iname (one entry per line)
// and store it in k2 format
// the compressed matrix is stored to :a and its size to :msize
// if :xsize>0 that value is forced to be the size of k2 matrix
// return the internal size of the k2 matrix (which has the form 2**k*MMsize)
// since entries are encoded in 64 bits, each index can be at most 32 bits
// so the maximum matrix size is 2^32 (change ia[],imin,imax type to go further)
size_t mread_from_textfile(size_t *msize, k2mat_t *a, char *iname, size_t xsize)
{
  assert(iname!=NULL && a!=NULL);
  FILE *f = fopen(iname,"rt");
  if(f==NULL) quit("mread_from_file: cannot open input file",__LINE__,__FILE__);
  // generate interleaved array from input file
  size_t n; // number of entries
  // since we are storing entries in 64 bits each index must fit in 32 bits   
  if(xsize>1UL+UINT32_MAX) quit("mread_from_textfile: matrix too large, current limit is 2^32",__LINE__,__FILE__);
  uint64_t *ia = create_ia(f,&n,msize,xsize);
  assert(xsize==0 || *msize==xsize);
  fclose(f);
  // compress the matrix represented by the ia[] array into a k2mat_t structure
  size_t asize = mread_from_ia(ia,n,*msize,a);
  free(ia);
  return asize;
}



// write the content of the :msize x :msize k2 matrix :a to a
// text file in one entry per line format
void mwrite_to_textfile(size_t msize, size_t asize, const k2mat_t *a, char *outname)
{
  assert(outname!=NULL && a!=NULL);
  assert(asize>=msize);
  FILE *f = fopen(outname,"wt");
  if(f==NULL) quit("mwrite_to_file: cannot open output file",__LINE__,__FILE__);

  if(k2is_empty(a)) {  // an empty k2 matrix has no entries
    fclose(f);
    return;
  }
  size_t pos = 0;
  mdecode_to_textfile(f,msize,0,0,asize,a,&pos);
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
// since entries are encoded in 64 bits, each index can be at most 32 bits
// so the maximum matrix size is 2^32 (change ia[] type to go further)
static size_t mread_from_ia(uint64_t ia[], size_t n, size_t msize, k2mat_t *a)
{
  assert(ia!=NULL && a!=NULL);
  assert(n>0);                // we cannot represent an empty matrix
  assert(msize>1);
  assert(a_eq_b2(n,msize) || a_lt_b2(n,msize));   // entries can be at most msize**2
  k2_free(a);                 // free previous content of a 
  if(msize>1UL+UINT32_MAX) quit("mread_from_ia: matrix too large, current limit is 2^32",__LINE__,__FILE__);
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
  // encode ia[0,n-1] into the k2mat_t structure a
  mencode_ia(ia,n,0,asize,a);
  return asize;
}

// compare a and b^2 with only operations
// involving uint64_t and without overflow 
static __inline__ bool a_eq_b2(uint64_t a, uint64_t b)
{
  return (a/b==b) ? (a%b==0) : false;
}

static __inline__ bool a_lt_b2(uint64_t a, uint64_t b) 
{
  return (a/b<b);
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
//   imin  smallest value assigned to the current submatrix 
//   size  submatrix size (has the form 2^k*MMsize)
//   *c    output k2mat_t structure to be filled in dfs order
// all entries in ia[0,n-1] are in the range [imin, imin+size*size) 
// all these entries must be encoded in the k2mat c 
// In previous versions of the code also the parameter imax = imin+size^2
// was used explicitly: it has been removed since for size==2^32
// such value could be 2^64 and tehrefore not representable in a uint64  
// called by mread_from_ia()
static void mencode_ia(uint64_t *ia, size_t n, uint64_t imin, size_t size, k2mat_t *c) {
  //printf("Size=%zu, n=%zu, imin=%lu\n",size,n,imin);
  assert(ia!=NULL);
  assert(n>0);
  assert(ia[0]>=imin); 
  // assert(ia[n-1]<imin+size*size); replaced by the following line
  assert( a_lt_b2(ia[n-1]-imin, size)); 
  assert(size%2==0 && size>=2*MMsize);
  // case of a full submatrix
  if(a_eq_b2(n,size)) { // equivalent to (n==size*size) but no overflow   
    k2add_node(c,ALL_ONES);       // submatrix is full 
    return;
  }
  // determine range of submatrices
  assert(size/2<UINT32_MAX);  // check that size/2 can be squared without overflow 
  uint64_t range = (size/2)*(size/2);
  uint64_t left = imin + range;
  uint64_t mid = left+range;
  uint64_t right = mid+range;
  // printf("range=%lu imax-imin=%lu\n",range,right+range-imin);
  if(size==1UL+UINT32_MAX) // max value size=2^32 treated separately
    assert(right-imin>0 && right-imin+range==0); // equiv to right-imin+range==2^64
  else   
    assert(a_eq_b2(right-imin+range,size));  // equiv to: right+range == imin + size^2
  // determine range in ia[] of the 4 submatrices entries
  size_t imid = binsearch(ia,n,mid);    //   first entry of A[10]
  size_t ileft = imid>0 ? binsearch(ia,imid,left):0; // first entry of A[01]
  size_t iright = imid<n? binsearch(ia+imid,n-imid,right)+imid:n; // first entry of A[11]
  // the four submatrices are: 
  //    ia[0,ileft-1], ia[ileft,imid-1], ia[imid,iright-1], ia[iright,n-1]
  // and contain values in the ranges
  //    [imin,left), [left,mid), [mid,right), [right,imin+size^2)
  // start building c
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_ONES as root placeholder 
  node_t rootc=NO_CHILDREN;                 // actual root node to be computed
  // here we are assuming that the submatrices are in the order 00,01,10,11

  if(ileft>0) { // submatrix 00 is not empty
    rootc |= (1<<0); // set 00 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia,ileft,imin,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia,ileft,imin,size/2);
      k2add_minimat(c,cx);
    }
  }

  if(ileft<imid) { // submatrix 01 is not empty
    rootc |= (1<<1); // set 01 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia+ileft,imid-ileft,left,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia+ileft,imid-ileft,left,size/2);
      k2add_minimat(c,cx);
    }
  }

  if(iright>imid) { // submatrix 10 is not empty
    rootc |= (1<<2); // set 10 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia+imid,iright-imid,mid,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia+imid,iright-imid,mid,size/2);
      k2add_minimat(c,cx);
    }
  }

  if(iright<n) { // submatrix 11 is not empty
    rootc |= (1<<3); // set 11 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia+iright,n-iright,right,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia+iright,n-iright,right,size/2);
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
// since entries are encoded in 64 bits, each index can be at most 32 bits
// so the maximum matrix size is 2^32 (change ia[] type to go further)
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize, size_t xsize)
{
  int64_t maxentry = 0; // largest entry in the file
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
    // since we are storing entries in 64 bits each index must fit in 32 bits       
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
  // save output parameters   
  if(xsize==0) { // if xsize==0 size is largest index + 1
    if(maxentry+1>SIZE_MAX)  // highly unlikely, but you never know... 
      quit("create_ia: cannot represent matrix size",__LINE__,__FILE__);
    *msize = (size_t) maxentry+1;
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
  assert(i<msize+2*size && j<msize+2*size);
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
  for(size_t k=0;k<4;k++) {  
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

