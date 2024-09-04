/* Routines for converting boolean matrices in text form
   to/from compressed k^2 tree representation 

   Note: internally matrix dimensions are always of the form 2^k times the size 
   of a minimatrix (those stored at the leaves of the tree), with k>0
   (somewhere this is called the k2_internal_size); the input can be of 
   any size (not larger than that) and the k2 matrix is padded with 0's 
   (virtually since they are not stored)

   The conversion txt->k2 is done using an auxiliary "interleaved" array:
   each matrix entry consists of two uint32_t (row and column indices).
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
   a_lt_b2 and a_eq_b2 testing whether a<b*b or a==b*b without multiplications 
    
   The conversion k2->txt is done doing a visit of the tree in preorder and
   each time a nonzero entry is found its indices are written to the output file 
   
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
static inline bool a_eq_b2(uint64_t a, uint64_t b);
static inline bool a_lt_b2(uint64_t a, uint64_t b);

// limit for saving subtree size
extern int32_t Depth_subtree_size_save;

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
// subtree size computation and verification 
int32_t Depth_subtree_size_save=0;


// resizable vector of uint64's
typedef struct {
  uint64_t *v;       // array of u64's
  size_t n;          // number of elements 
  size_t nmax;       // maximumt capacity 
} vu64_t;


void vu64_init(vu64_t *z)
{
  z->nmax=16;
  z->n=0;
  z->v = malloc(z->nmax*sizeof *(z->v) );
}

void vu64_free(vu64_t *z)
{
  free(z->v);
  z->v = NULL; 
  z->nmax = z->n=0;
}

// add i elements at the end of z
void vu64_grow(vu64_t *z, size_t i) 
{
  z->n +=i;
  if(z->n>z->nmax) {
    z->nmax *= 2;
    z->v = realloc(z->v, z->nmax*sizeof *(z->v) );
    if(z->v==NULL) quit("realloc failed",__LINE__,__FILE__);
  }
}


// The next function does a dfs visit of the k2 matrix :m and write 
// the subtree sizes and the econding of the sizes in the growing array z
// the returned value is the total size of the k2 submatrix (in nibbles) 
// and in the upper 24 bits the cost of subtree encoding (ie the cost
// of encoding recursively the submatrices) 
// 
// Explanation:
// Assuming the root of T has tree children the subtree represention is something like:
//   R111111111222222333333333333
// we need to compute and return the total size of the representation ie the length of 
// the above string (ie #R(==1) + #1 + #2 + #3)  and we need to store in z an encoding 
// of #1 and #2 followed by the same information for the subtrees 1, 2 and 3
// This information stored in z for T will be called the subtree information for T.
// We can fill z with a DFS visit of the k2tree, when we reach the a subtree T
// aboe we leave two empty slots in z, call the function recursively getting 
// #1, #2 and #3, storing #1, #2 in the empty slots and return 1 + #1 + #2 + #3
// However, since z is used to skip subtree 1 and/or 2, in z together with 
// #1 we also need to store the total information stored in z for the subtree 
// rooted at 1, and the same for the subtree rooted at 2.
// Hence, array t will not simply contain the encoding of 
//     <#1> <#2> info_Sub(1) info_Sub(2) info_Sub(3)
// (where < > denotes an encoding, for example 7x8), but
//     <#1> <|info_Sub(1)|> <#2> <|info_Sub(2)|> info_Sub(1) info_Sub(2) info_Sub(3) 
// (to complicate things info Sub(i) is different from above because 
//  now includes the additional values <|info Sub()|> for the subtrees  
// To do the computation, this function, for each subtree T, returns 2 values:
//   1. the total length of T encoding (as before, 1+ #1 + #2 + #3)
//   2. the total length of the above complete encoding of T subtree information
//      If E1, E2, E3 are the lenghts of the encodings 
//      for the subtree  Ei = |info Sub(i)|   (obtained by the recursive calls) 
//      then the cost of the encoding of T is 
//        |<#1>| + |<E1>| + |<#2>| + |<E2>| + E1 + E2 + E3
//      note that in z's empty slots the function has to store 
//      the values #1 E1 #2 E2
// To simplify the code (!!!) the function returns a single uint64, with
// the less significant 40 bits storing the total length of T (item 1)
// and in the more significant 24 bits the lengths of the econdings
// Hence the recursive calls return: (here / means justapoxition) 
//      A1 = #1/E1,    A2 = #2/E2,    A3 = #3/E3
// and the function should return
//      1+ #1 + #2 + #3 / |<#1>| + |<E1>| + |<#2>| + |<|E2|>| + E1 + E2 + E3
// assuming there are no overflows, the desired value is
//      1+A1+A2+A3 + (|<#1>| + |<E1>| + |<#2>| + |<|E2|>|)<<48
// The above scheme is valid for ordinary internal nodes. 
// If T is an ALL_ONES leaf, then its size is 1 and the size of the 
//  subtree sizes encoding is 0.
// If T has height 1, then its size is 1+#child*Minimat_node_ratio
//  but there is no need to store the subtree size information since 
//  each subtree has size Minimat_node_ratio
// If T has depth2go==1 we need to store the subtree size information 
//  for T as usual, but we know that we are not storing information for
//  T subtrees so E1=E2=E3=0. Since this is something we can check
//  during the visit we can simply avoid storing <E1> and <E2>
//  (at the moment for simplicity we do store them)
// If T has depth2go<=0 we need to report T size as usual, but 
//  we do not need to store any subtree information and we report 0 
//  as the total lenght of the subtree information
// Note that in this function we are not actually encoding the values 
// but computing the values that will be later encoded. Such values are stored 
// in z using the above 40+24 scheme, the values then have to be stored (on disk)
// using the appropriate scheme. In a first attempt we can avoid the encoding
// and just use the array z as above. In that case we can measure everything
// in uint64's so we simply have that |<#1>| + |<E1>| = 1 (1 uint64_t)
//
// Note: the above scheme can be probably improved in speed with a minimal
// space increase. Given the structure 
//  <#1> <|info_Sub(1)|> <#2> <|info_Sub(2)|> info_Sub(1) info_Sub(2) info_Sub(3) 
// a major issue is that to reach the info_Sub() information one has to
// to skip the "<#1> <|info_Sub(1)|> <#2> <|info_Sub(2)|>" part and then
// possibly sum together |info_Sub(1)|, |info_Sub(2)| etc. So we could
// store also len(<#1> <|info_Sub(1)|> <#2> <|info_Sub(2)|>), and we could 
// store |info_Sub(1)| + |info_Sub(2)| instead of |info_Sub(2)| and so on. 

// Constant to store size and esizes in a single value
#define BITSxTSIZE 40
#define TSIZEMASK ( (((uint64_t) 1)<<BITSxTSIZE) -1 )
// note: potential overflow if tree sizes cannot be expressed in BITSxTSIZE bits
// and if the encoding size cannot be expressed with 64-BITSxTSIZE 
// setting BITSxTSIZE at least 40 make the first event unlikely,
// while the second event is possible if we keep information for many levels. 
// The first event should be detected by the test on *pos immediately 
// before the final return. The second event is detected using 
// __builtin_add_overflow (-ftrapv or fsanitize do not work since they are 
// for signed int and they would add extra checks for all operations)
// 
#define CHECK_ESIZE_OVERFLOW 1
// Parameters:
//  size of the current submatrix
//  m,*pos the current submatrix starts at position *pos within *m 
//  z dynamic vector where the subtree information will be stored
//  # levels for which we store the subtree information 
uint64_t k2dfs_sizes(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, int32_t depth2go)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos); // implies m is non-empty
  size_t pos_save = *pos;  // save starting position of subtree
  node_t root = k2read_node(m,*pos); (*pos)++;
  assert(root<ILLEGAL_NODE); 
  if(root==ALL_ONES) { 
    assert(Use_all_ones_node); // all 1's matrix consists of root only, 
    return 1;                  // size is 1 subtree encoding size 0
  }
  // we have a non-singleton subtree to traverse 
  // compute number of children
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);

  size_t zn_save = z->n; // save starting position in size_array[]
  if(depth2go>0)
    vu64_grow(z,nchildren-1);
  size_t subtree_size = 1;     // account for root node
  size_t child_size = 0;       // size/esize of a child subtrees
  size_t csize[4];             // sizes/esizes of the children subtrees
  size_t cpos = 0;             // current position in size[]
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      if(size==2*MMsize)  // end of recursion
        *pos += (child_size = Minimat_node_ratio);
      else { // recurse on submatrix
        child_size =  k2dfs_sizes(size/2,m,pos,z,depth2go-1); // read submatrix and advance pos
      }
      #ifdef CHECK_ESIZE_OVERFLOW
      // save size and esize for possible later storage in z
      csize[cpos++] = child_size;
      // sum sizes and esizes, check for possible overflow    
      if(__builtin_add_overflow(subtree_size,child_size,&subtree_size))
        quit("Overflow in subtree encoding: make BITSxTSIZE smaller if possible",__LINE__,__FILE__);      
      #else
      subtree_size += (csize[cpos++] = child_size); // save and sum sizes and esizes, no check
      #endif
    }
  assert(cpos==nchildren); // we should have visited all children
  // check subtree size for all children except last one
  if(depth2go>0) {
    for(int i=0; i<cpos-1; i++)
      z->v[zn_save++] = csize[i];
    #ifdef CHECK_ESIZE_OVERFLOW
      if(__builtin_add_overflow(subtree_size,(nchildren-1)<<BITSxTSIZE,&subtree_size))
        quit("Overflow in subtree encoding: make BITSxTSIZE smaller if possible",__LINE__,__FILE__);      
    #else
    subtree_size += (nchildren-1)<<BITSxTSIZE;
    #endif
  }
  else assert(subtree_size>>BITSxTSIZE == 0); // there should not be any subtree encodings
  if(*pos != pos_save + (subtree_size&TSIZEMASK)) { // double check size
    fprintf(stderr,"Scanned size: %zu, computed size: %lu\n", *pos-pos_save,subtree_size&TSIZEMASK); 
    quit("Error or overflow in size encoding",__LINE__,__FILE__);
  } 
  return subtree_size;
}


// do a dfs visit of the k2 matrix :m and make sure subtree sizes match the ones in z
// the checking is done recursively up to depth depth2check 
static size_t k2dfs_check_sizes_depth(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, int32_t depth2check)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos);        // implies m is non-empty
  size_t pos_save = *pos;     // save starting position of subtree
  node_t root = k2read_node(m,*pos); (*pos)++;
  assert(root<ILLEGAL_NODE);
  if(root==ALL_ONES) { 
    assert(Use_all_ones_node); // all 1's matrix consists of root only, 
    return 1;                  // subtree size is 1 no subtree info to check
  }
  // we have a non-singleton subtree T to traverse 
  // compute number of children
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);
  // read subtree information if available 
  uint64_t *subtree_info = NULL;
  if(depth2check>0) {
    subtree_info = &(z->v[z->n]);  // starting point of T information
    z->n += nchildren-1;
  }
  size_t cnum = 0;             // current position in size[]
  size_t subtree_size = 1;     // account for root node
  for(int i=0;i<4;i++) 
    // invariant: both *pos and z->n point to the beginning of subtree cnum
    if(root & (1<<i)) {
      size_t pc = *pos;
      size_t nc = z->n;
      size_t child_size;
      if(size==2*MMsize)  // end of recursion
        *pos += (child_size = Minimat_node_ratio);
      else { // recurse on submatrix
        child_size =  k2dfs_check_sizes_depth(size/2,m,pos,z,depth2check-1);
      }
      assert(pc + child_size == *pos);  // check *pos
      if(depth2check>0 && cnum<nchildren-1) { // chek correctness of info in z
        assert((subtree_info[cnum]&TSIZEMASK) == child_size);
        // assert((subtree_info[cnum]>>BITSxTSIZE) == z->n - nc);
        if((subtree_info[cnum]>>BITSxTSIZE) != z->n - nc) 
          printf("Size: %zu, cnum %zu info:%ld, deltan: %zu\n",size,cnum, subtree_info[cnum]>>BITSxTSIZE, z->n-nc);
      }
      cnum++;
      subtree_size += child_size;
    }  
  assert(cnum==nchildren); // we should have visited all children
  assert(*pos == pos_save + subtree_size); // check subtree size
  return subtree_size;
}


// do a dfs visit of the k2 matrix :m and make sure subtree sizes match the ones in z
// the checking is done recursively until the subtree encodings are >0
// here we are calling tree the one we are exploring, 
// and subtrees its immediate descendant
// if check is true the the encoding of the current tree is stored in z
// and has to be checked. Recall that if the tree has 3 non empty children
// its encoding consists of 
//  <T1> <Sub1> <T2><Sub2> Sub1 Sub2 Sub3 (where <x> denotes size of x) 
static size_t k2dfs_check_sizes(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, 
                                size_t tot_encode_size)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos);        // implies m is non-empty
  size_t pos_save = *pos;     // save starting position of tree
  node_t root = k2read_node(m,*pos); (*pos)++;
  assert(root<ILLEGAL_NODE);
  if(root==ALL_ONES) { 
    assert(Use_all_ones_node); // all 1's matrix consists of root only, 
    return 1;                  // tree size is 1 no subtree info to check
  }
  // we have a non-singleton tree T to traverse 
  // compute number of children
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);
  // read subtree information if available 
  uint64_t *subtree_info = &(z->v[z->n]); // array with subtree_info information
  z->n += nchildren-1;                    // advance z->n to subtree encoding
  size_t encode_seen = nchildren-1;       // consume one item x non-last children 
  // visit children 
  size_t cnum = 0;             // current child 
  size_t tree_size = 1;        // account for root node
  for(int i=0;i<4;i++) 
    // invariant: both *pos and z->n point to the beginning of subtree cnum
    if(root & (1<<i)) {
      size_t pc = *pos;  // save current position in m and z 
      size_t nc = z->n;
      size_t child_subtree_size=0;
      size_t child_encode_size = 0;
      // compute size of subtree encoding
      if(cnum<nchildren-1) {
        child_encode_size = subtree_info[cnum]>>BITSxTSIZE;
        encode_seen += subtree_info[cnum]>>BITSxTSIZE;
      }
      else  // last child 
        child_encode_size = tot_encode_size -encode_seen; // remaining encoding 
      // ------- go down one level ----------- 
      if(size==2*MMsize)  // end of recursion
        *pos += (child_subtree_size = Minimat_node_ratio);
      else if(child_encode_size==0) 
        k2dfs_visit_fast(size/2,m,pos);  // advance pos to the of subtree 
      else // recurse on subtree
        child_subtree_size =  k2dfs_check_sizes(size/2,m,pos,z,child_encode_size);
      // check that child_subtree_size matches the advancement in *pos
      if(child_subtree_size != *pos -pc) 
        fprintf(stderr,"Subtree scanned size: %zu, reported size: %zu\n",*pos-pc,child_subtree_size);
      // if not last child check that stored subtree size matches
      if(cnum<nchildren-1 && child_subtree_size!=(subtree_info[cnum]&TSIZEMASK))
        fprintf(stderr,"Subtree stored size: %zu, reported size: %zu\n",subtree_info[cnum]&TSIZEMASK,child_subtree_size);
      // check stored subtree encoding matches 
      size_t scanned_encoding = z->n-nc;
      if(child_encode_size!=scanned_encoding) {
        if(cnum<nchildren-1) 
          fprintf(stderr,"Subtree encoding stored size: %zu, scanned size: %zu\n",child_encode_size,scanned_encoding);
        else
          fprintf(stderr,"Subtree encoding computed size: %zu, scanned size: %zu\n",child_encode_size,scanned_encoding);
      }
      cnum++;
      tree_size += child_subtree_size;
    }  
  assert(cnum==nchildren); // we should have visited all children
  assert(*pos == pos_save + tree_size); // check again tree size
  return tree_size;
}




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
  if(Depth_subtree_size_save > 0) {
    puts("Computing and checking subtree sizes");
    vu64_t z;
    vu64_init(&z);
    size_t pos=0,znsave;
    uint64_t p = k2dfs_sizes(asize,a,&pos,&z,Depth_subtree_size_save);
    assert(pos==a->pos);
    printf("Returned by k2dfs_sizes: %zu, sizes stored: %zu\n", p&TSIZEMASK, z.n);
    pos=0;
    znsave = z.n;
    z.n=0;
    size_t pcheck = k2dfs_check_sizes(asize,a,&pos,&z,znsave);
    printf("Returned by k2dfs_check_sizes: %zu, sizes read: %zu\n", pcheck, z.n);
    vu64_free(&z);
  }
  return asize;
}

// compare a and b^2 with only operations
// involving uint64_t and without overflow 
static inline bool a_eq_b2(uint64_t a, uint64_t b)
{
  return (a/b==b) ? (a%b==0) : false;
}

static inline bool a_lt_b2(uint64_t a, uint64_t b) 
{
  return (a/b<b);
}

// given a sorted uint64_t array ia[0,n-1] containing distinct values find 
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
//   smin  smallest value assigned to the current submatrix 
//   size  submatrix size (has the form 2^k*MMsize)
//   *c    output k2mat_t structure to be filled in dfs order
// all entries in ia[0,n-1] are in the range [smin, smin+size*size) 
// all these entries must be encoded in the k2mat c 
// In previous versions of the code also the parameter imax = smin+size^2
// was used explicitly: it has been removed since for size==2^32
// such value could be 2^64 and therefore not representable in a uint64  
// called by mread_from_ia()
static void mencode_ia(uint64_t *ia, size_t n, uint64_t smin, size_t size, k2mat_t *c) {
  //printf("Size=%zu, n=%zu, smin=%lu\n",size,n,smin);
  assert(ia!=NULL);
  assert(n>0);
  assert(ia[0]>=smin); 
  // assert(ia[n-1]<smin+size*size); replaced by the following line
  assert( a_lt_b2(ia[n-1]-smin, size)); 
  assert(size%2==0 && size>=2*MMsize);
  // case of a full submatrix
  if(a_eq_b2(n,size) && Use_all_ones_node) { // equivalent to (n==size*size) but no overflow   
    k2add_node(c,ALL_ONES);       // submatrix is full 
    return;
  }
  // determine range of submatrices
  assert(size/2<UINT32_MAX);  // check that size/2 can be squared without overflow 
  uint64_t range = (size/2)*(size/2);
  uint64_t left = smin + range;
  uint64_t mid = left+range;
  uint64_t right = mid+range;
  // printf("range=%lu imax-smin=%lu\n",range,right+range-smin);
  if(size==1UL+UINT32_MAX) // max value size=2^32 treated separately
    assert(right-smin>0 && right-smin+range==0); // equiv to right-smin+range==2^64
  else   
    assert(a_eq_b2(right-smin+range,size));  // equiv to: right+range == smin + size^2
  // determine range in ia[] of the 4 submatrices entries
  size_t imid = binsearch(ia,n,mid);    //   first entry of A[10]
  size_t ileft = imid>0 ? binsearch(ia,imid,left):0; // first entry of A[01]
  size_t iright = imid<n? binsearch(ia+imid,n-imid,right)+imid:n; // first entry of A[11]
  // the four submatrices are: 
  //    ia[0,ileft-1], ia[ileft,imid-1], ia[imid,iright-1], ia[iright,n-1]
  // and contain values in the ranges
  //    [smin,left), [left,mid), [mid,right), [right,smin+size^2)
  // start building c
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_ONES as root placeholder 
  node_t rootc=NO_CHILDREN;                 // actual root node to be computed
  // here we are assuming that the submatrices are in the order 00,01,10,11

  if(ileft>0) { // submatrix 00 is not empty
    rootc |= (1<<0); // set 00 bit
    if(size>2*MMsize) // if size>2*MMsize recurse
      mencode_ia(ia,ileft,smin,size/2,c);
    else { // size==2*MMsize: write a minimatrix
      minimat_t cx = minimat_from_ia(ia,ileft,smin,size/2);
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


