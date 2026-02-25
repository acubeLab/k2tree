/* Routines for converting boolean matrices in sparse text form
   (ie one entry per line) to/from compressed k^2 tree representation 

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
#include "minimats.c" // includes k2.h bbm.h

#include "libsais/include/libsais64.h"
#include "pointers.h"
#include "rank_0000.h"


// prototypes of static functions
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize, size_t xsize);
static size_t mread_from_ia(uint64_t ia[], size_t n, size_t msize, k2mat_t *a);
static void mencode_ia(uint64_t *ia, size_t n, uint64_t imin, size_t size, k2mat_t *c);
static void mdecode_to_textfile(FILE *outfile, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos);
static size_t binsearch(uint64_t *ia, size_t n, uint64_t x);
static inline bool a_eq_b2(uint64_t a, uint64_t b);
static inline bool a_lt_b2(uint64_t a, uint64_t b);


// read a matrix from the text file :iname (one entry per line)
// and store it in k2 format
// the compressed matrix is stored to :a and its size to a->realsize
// if :xsize>0 that value is forced to be the size of k2 matrix
// return the actual size of the k2 matrix
// since entries are encoded in 64 bits, each index can be at most 32 bits
// so the maximum matrix size is 2^32 (change ia[],imin,imax type to go further)
size_t mread_from_textfile(k2mat_t *a, char *iname, size_t xsize)
{
  assert(iname!=NULL && a!=NULL);
  FILE *f = fopen(iname,"rt");
  if(f==NULL) quit("mread_from_textfile: cannot open input file",__LINE__,__FILE__);
  // generate interleaved array from input file
  size_t n; // number of entries
  size_t msize; // computed matrix size;
  // since we are storing entries in 64 bits each index must fit in 32 bits   
  if(xsize>1UL+UINT32_MAX) quit("mread_from_textfile: matrix too large, current limit is 2^32",__LINE__,__FILE__);
  uint64_t *ia = create_ia(f,&n,&msize,xsize);
  assert(xsize==0 || msize==xsize);
  fclose(f);
  // compress the matrix represented by the ia[] array into a k2mat_t structure
  size_t asize = mread_from_ia(ia,n,msize,a);
  free(ia);
  a->realsize = msize; a->fullsize = asize;
  return msize; // return the real size of the k2_matrix
}



// write the content of the :msize x :msize k2 matrix :a to a
// text file in one entry per line format
// respcet main_diag_1 and backp pointers. subtinfo is not used
void mwrite_to_textfile(const k2mat_t *a, char *outname)
{
  assert(outname!=NULL && a!=NULL);
  size_t msize = a->realsize, asize = a->fullsize; 
  assert(asize>=msize);
  FILE *f = fopen(outname,"wt");
  if(f==NULL) quit("mwrite_to_file: cannot open output file",__LINE__,__FILE__);

  if(k2is_empty(a)) {  // an empty k2 matrix has no entries
    if(a->main_diag_1) { // but if main_diag_1 is on we have the diagonal entries
      for(size_t i=0;i<msize;i++) {
        int e = fprintf(f,"%zu %zu\n",i,i);
        if(e<0) quit("mwrite_to_textfile: error writing to output file",__LINE__,__FILE__);
      }
    }
    fclose(f);
    return;
  }
  size_t pos = 0;
  mdecode_to_textfile(f,msize,0,0,asize,a,&pos);
  fclose(f);
  assert(pos==k2pos(a)); // check we read all the k2mat_t structure
}

// subtree size computation and verification 

// The next function does a dfs visit of the k2 matrix :m and writes 
// the subtree sizes and the enconding of the sizes in the growing array z
// The returned value is the total size of the k2 submatrix (in nibbles) 
// and in the upper 24 bits the cost of subtree encoding (ie the cost
// of encoding recursively the submatrices) 
// 
// Explanation:
// Assuming the root R of T has three children the subtree representation is something like:
//   R111111111222222333333333333
// we need to compute and return the total size of the representation ie the length of 
// the above string (ie #R(==1) + #1 + #2 + #3)  and we need to store to z an encoding 
// of #1 and #2 followed by the same information for the subtrees 1, 2 and 3
// This information stored in z is called the subtree information for T.
// We can fill z with a DFS visit of the k2tree: when we reach a subtree T
// above we leave two empty slots in z, call the function recursively getting 
// #1, #2 and #3, store #1, #2 in the empty slots and return 1 + #1 + #2 + #3
// However, since z is used to skip subtree 1 and/or 2, in z together with 
// #1 we also need to store the total information stored in z for the subtree 
// rooted at 1, and the same for the subtree rooted at 2 (because we need
// to be able to skip this information as well). 
// Hence, array z will not simply contain the encoding of 
//     <#1> <#2> info_Sub(1) info_Sub(2) info_Sub(3)
// (where < > denotes an encoding of a size, for example 7x8 byte encoding), but
//     <#1> <|info_Sub(1)|> <#2> <|info_Sub(2)|> info_Sub(1) info_Sub(2) info_Sub(3) 
// (to complicate things info Sub(i) is different from above because 
//  now includes the additional values <|info Sub()|> for the subtrees of i  
// To do the computation for each subtree T  this function returns 2 values:
//   1. the total length of T encoding (as before, 1+ #1 + #2 + #3)
//      this is the amount of info that we need to skip in the k2 matrix to skip T
//   2. the total length of the above complete encoding of T subtree information
//      this is the amont of info we need to skip in z to skip T 
//      If E1, E2, E3 are the lenghts of the encodings 
//      for the subtree  Ei = |info Sub(i)|   (obtained by the recursive calls)
//      in z for T we store  
//         <#1> <E1> <#2> <E2> info_Sub(1) info_Sub(2) info_Sub(3)
//      hence the cost of the encoding (of the subtrees) of T is 
//        |<#1>| + |<E1>| + |<#2>| + |<E2>| + E1 + E2 + E3
//      note that in the above example in z's empty slots the function has 
//      to store the values #1 E1, #2 E2
// To simplify the code (!!!) the function returns a single uint64, with
// the less significant 40 bits storing the total length of T (item 1)
// and in the more significant 24 bits the lengths of the encodings
// Hence the recursive calls return: (here / means justapoxition lowbits/hibits) 
//      A1 = #1/E1,    A2 = #2/E2,    A3 = #3/E3
// and the function should return
//      1+ #1 + #2 + #3 / |<#1>| + |<E1>| + |<#2>| + |<|E2|>| + E1 + E2 + E3
// assuming there are no overflows, the desired value is
//      1+A1+A2+A3 + (|<#1>| + |<E1>| + |<#2>| + |<|E2|>|)<<40
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
//  (at the moment we do store them, because it allows us
//   to use more complex schemes, see below "An alternative...")
// If T has depth2go<=0 we need to report T size as usual, but 
//  we do not need to store any subtree information and we report 0 
//  as the total lenght of the subtree information
// 
// An alternative scheme is to store the subtree info only for 
// those subtrees which are larger than a certain threshold;
// in this case we need a strategy to recognize, when we use the matrix, if
// the subtree info is present without additional external information.
// The idea is the following (see k2split_k2() for an example):
//   if m->subtinfo==NULL then there is no subtree info for the current 
//      tree (and its subtrees!) 
//   if m->subtinfo!=NULL then it contains the size of its subtrees,
//    in the above example #1 and #2, the size of the last subtree 
//    is obtained by subtraction: #3 = #T -1 (root) - #1 -#2
//   for each subtree (1,2,3) we need to compute Ei to see if there
//    is subtree info stored for Ei (ie if Ei!=0)
//    for 1 and 2 we just look at <E1>, <E2> that are saved in subtinfo
//    while we get E3 by the formula 
//       m->subtinfo_size = |<#1>| + |<E1>| + |<#2>| + |<E2>| + E1 + E2 + E3
//    with our simple encoding it is |<#1>| + |<E1>| + |<#2>| + |<E2>| = 2
//    hence:
//       E3 = m->subtinfo_size - (nchildren-1) - E1 - E2
// Note that this requires that each time we create a (sub)matrix we 
// maintain the correct m->subtinfo_size (which has no other uses).

// Note that in this function we are not actually encoding (ie compressing) the values 
// but only computing them (we can always encode later). Such values are stored 
// to z using the above 40+24 scheme, the values then have to be stored (on disk)
// using the appropriate scheme. As a first attempt we avoid the encoding
// and just use the array z as above. In that case we can measure everything
// in uint64's so we simply have that |<#1>| + |<E1>| = 1 (1 uint64_t)
//
// Note that if we are going to use a more complex encoding (say 7x8)
// we cannot easily fill the array z from left to right since the 
// the size of the empty slots we create at the beginning of T's visit
// is unknow at the beginnig of the visit since it depends on the 
// subtrees content. A possibile solution could be a 2-pass encoding
// (in the first pass we compute the correct sizes but we store
// them in uint64s, then we do the actual encoding using say 7x8,
// the drawback is larger working space), or we fill z in reverse
// (first the subtrees and then T, visiting the tree right to left)
// and then we write it to disk in reverse, drawback: complex code)   

// Note: the above scheme can be probably improved in speed with a minimal
// space increase. Given the structure 
//  <#1> <|info_Sub(1)|> <#2> <|info_Sub(2)|> info_Sub(1) info_Sub(2) info_Sub(3) 
// a major issue is that to reach the info_Sub() information one has to
// to skip the "<#1> <|info_Sub(1)|> <#2> <|info_Sub(2)|>" part and then
// possibly sum together |info_Sub(1)|, |info_Sub(2)| etc. So we could
// store also len(<#1> <|info_Sub(1)|> <#2> <|info_Sub(2)|>), and we could 
// store |info_Sub(1)| + |info_Sub(2)| instead of |info_Sub(2)| and so on. 

// Constants to store size and esizes in a single value (moved to k2.h) 
// #define BITSxTSIZE 40
// #define TSIZEMASK ( (((uint64_t) 1)<<BITSxTSIZE) -1 )
// note: potential overflow if tree sizes cannot be expressed in BITSxTSIZE bits
// and if the encoding size cannot be expressed with 64-BITSxTSIZE 
// setting BITSxTSIZE at least 40 make the first event unlikely,
// while the second event is possible if we keep information for many levels. 
// The first event should be detected by the test on *pos immediately 
// before the final return. The second event is detected using 
// __builtin_add_overflow (-ftrapv or -fsanitize do not work since they are 
// for signed int and they would add extra checks for all operations)
// 
#define CHECK_ESIZE_OVERFLOW 1
// The next function does a dfs visit of the k2 matrix :m and writes 
// the subtree sizes and the enconding of the sizes in the growing array z
// The returned value is the total size of the k2 submatrix (in nibbles) 
// and in the upper 24 bits the cost of subtree encoding (ie the cost
// of encoding recursively the submatrices) 
// 
// Parameters:
//  size: k2 size of the current submatrix (ie 2^k * MMsize)
//  m,*pos the current submatrix starts at position *pos within *m 
//  z: dynamic vector where the subtree information will be stored
//  depth2go: # levels for which we store the subtree information 
// uses the simpler .sinfo format with the last child  of each node excluded
// This function (with a depth limit) can be used only for the case in which 
// the subtreesize info for the last child is not stored (SIMPLEBACKPOINTERS case)
#ifdef SIMPLEBACKPOINTERS
uint64_t k2dfs_sizes(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, int32_t depth2go)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos); // implies m is non-empty
  size_t pos_save = *pos;  // save starting position of subtree
  node_t root = k2read_node(m,*pos); (*pos)++;
  assert(root<ILLEGAL_NODE); 
  if(root==ALL_ONES)            // all 1's matrix consists of root only, 
    return 1;                  // size is 1 subtree encoding size 0
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
  // add subtree size for all children except last one
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
  else assert(subtree_size>>BITSxTSIZE == 0); // there should not be any subtree encoding
  if(*pos != pos_save + (subtree_size&TSIZEMASK)) { // double check size
    fprintf(stderr,"Scanned size: %zu, computed size: %lu\n", *pos-pos_save,subtree_size&TSIZEMASK); 
    quit("Error or overflow in size encoding",__LINE__,__FILE__);
  } 
  return subtree_size;
}
#endif

// compute subtree information as above, but information is stored only 
// for large trees, ie when the number of nodes is larger than :limit
#ifdef SIMPLEBACKPOINTERS
uint64_t k2dfs_sizes_limit(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, size_t limit)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos); // implies m is non-empty
  size_t pos_save = *pos;  // save starting position of subtree
  node_t root = k2read_node(m,*pos); (*pos)++;
  assert(root<ILLEGAL_NODE); 
  if(root==ALL_ONES)           // all 1's matrix consists of root only, 
    return 1;                  // size is 1, subtree encoding size 0
  // we have a non-singleton subtree to traverse 
  // compute number of children
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);

  if(size==2*MMsize) { // end of recursion 
    *pos += nchildren*Minimat_node_ratio;
    return 1 +  nchildren*Minimat_node_ratio; // size, subtree encoding size 0 
  }

  size_t zn_save = z->n; // save starting position in size_array[]
  vu64_grow(z,nchildren-1);    // reserve space for children
  size_t subtree_size = 1;     // account for root node
  size_t child_size = 0;       // size/esize of a child subtrees
  size_t csize[4];             // sizes/esizes of the children subtrees
  size_t cpos = 0;             // current position in size[]
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      // recurse on submatrix
      child_size =  k2dfs_sizes_limit(size/2,m,pos,z,limit); // read submatrix and advance pos
      assert(child_size>0);
      // save size and esize for possible later storage in z
      csize[cpos++] = child_size;
      // add sizes and esizes to subtree_size, check for possible overflow    
      if(__builtin_add_overflow(subtree_size,child_size,&subtree_size))
        quit("Overflow in subtree encoding: make BITSxTSIZE smaller if possible",__LINE__,__FILE__);      
    }
  assert(cpos==nchildren); // we should have visited all children
  // check subtree size for all children except last one
  if((subtree_size&TSIZEMASK)>limit) {
    for(int i=0; i<cpos-1; i++)
      z->v[zn_save++] = csize[i];
    // add nchildren-1 to the encoding size, checking for overflow   
    if(__builtin_add_overflow(subtree_size,(nchildren-1)<<BITSxTSIZE,&subtree_size))
      quit("Overflow in subtree encoding: make BITSxTSIZE smaller if possible",__LINE__,__FILE__);      
  }
  else {
    assert(subtree_size>>BITSxTSIZE == 0); // there should not be any subtree encodings
    assert(z->n == zn_save+nchildren-1);   // no subtree information stored
    z->n = zn_save;                        // no subtree information
  }
  if(*pos != pos_save + (subtree_size&TSIZEMASK)) { // double check size
    fprintf(stderr,"Scanned size: %zu, computed size: %lu\n", *pos-pos_save,subtree_size&TSIZEMASK); 
    quit("Error or overflow in size encoding",__LINE__,__FILE__);
  } 
  return subtree_size;
}
#else
// NEW: version storing the info also for the last child  
uint64_t k2dfs_sizes_limit(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, size_t limit)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos); // implies m is non-empty
  size_t pos_save = *pos;  // save starting position of subtree
  node_t root = k2read_node(m,*pos); (*pos)++;
  assert(root<ILLEGAL_NODE); 
  if(root==ALL_ONES)           // all 1's matrix consists of root only, 
    return 1;                  // size is 1, subtree encoding size 0
  // we have a non-singleton subtree to traverse 
  // compute number of children
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);

  size_t zn_save = z->n; // save starting position in size_array[]
  vu64_grow(z,nchildren);
  size_t subtree_size = 1;     // account for root node
  size_t child_size = 0;       // size/esize of a child subtrees
  size_t csize[4];             // sizes/esizes of the children subtrees
  size_t cpos = 0;             // current position in size[]
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      if(size==2*MMsize)  // end of recursion
        *pos += (child_size = Minimat_node_ratio);
      else { // recurse on submatrix
        child_size =  k2dfs_sizes_limit(size/2,m,pos,z,limit); // read submatrix and advance pos
      }
      // save size and esize for possible later storage in z
      csize[cpos++] = child_size;
      // add sizes and esizes to subtree_size, check for possible overflow    
      if(__builtin_add_overflow(subtree_size,child_size,&subtree_size))
        quit("Overflow in subtree encoding: make BITSxTSIZE smaller if possible",__LINE__,__FILE__);      
    }
  assert(cpos==nchildren); // we should have visited all children
  // check subtree size for all children including the last one
  if((subtree_size&TSIZEMASK)>limit) {
    // complete encodig for current subtree with the x child info   
    for(int i=0; i<cpos; i++)
      z->v[zn_save++] = csize[i];
    // add nchildren to the encoding size, checking for overflow   
    if(__builtin_add_overflow(subtree_size,(nchildren)<<BITSxTSIZE,&subtree_size))
      quit("Overflow in subtree encoding: make BITSxTSIZE smaller if possible",__LINE__,__FILE__);      
  }
  else {
    assert(subtree_size>>BITSxTSIZE == 0); // there should not be any subtree encodings
    assert(z->n == zn_save+nchildren);     // no subtree information stored
    z->n = zn_save;                        // no subtree information
  }
  if(*pos != pos_save + (subtree_size&TSIZEMASK)) { // double check size
    fprintf(stderr,"Scanned size: %zu, computed size: %lu\n", *pos-pos_save,subtree_size&TSIZEMASK); 
    quit("Error or overflow in size encoding",__LINE__,__FILE__);
  } 
  return subtree_size;
}
#endif

// do a dfs visit of the k2 matrix :m and make sure subtree sizes match the ones in :z
// the checking is done recursively, but as soon as the encoding of a subtree
// has length 0, that subtree is explored with a fast dfs visit that only reports 
// the subtree size.
// In the code below we call "tree" the one we are exploring (representing :m) 
// and "subtrees" its immediate descendant
// Recall that if the tree has 3 non empty children
// its encoding consists of 
//   <T1> <Sub1> <T2> <Sub2> Sub1 Sub2 Sub3 (where <x> denotes size of x) 
// We compare <T1> and <T2> with the size returned from the subtree visits
// <T3> is not stored so it is checked at the upper level where
// 1 + <T1> + <T2> + <T3> will be compared with the size stored for 
// T's parent. <Sub1> and <Sub2> are tested with the amount of data
// scanned during the visit of T1 and T2. The value <Sub3> is obtained as
//  <Sub3> = tot_encode_size - <<T1> <Sub1> <T2><Sub2>> - <Sub1> - <Sub2>
// and is compared with the amount of data scanned during the visit of T3
// Parameters:
//  size  internal size of the current submatrix
//  m,*pos the current submatrix starts at position *pos within *m 
//  z  dynamic vector where the subtree information to be checked is stored
//  tot_encode_size total size of the encoding (info in z) for tree (and subtrees)
//                  as obtained at the previous level (T's parent)
//  Note: it is assumed that the root m[*pos] node has associate subtinfo in z[z->n]
// Return:
// size (number of nodes) of T (hence not including the info in z) 
#ifdef SIMPLEBACKPOINTERS
size_t k2dfs_check_sizes(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, 
                                size_t tot_encode_size)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos);        // implies m is non-empty
  size_t pos_save = *pos;     // save starting position of tree
  node_t root = k2read_node(m,*pos); (*pos)++;
  assert(root<ILLEGAL_NODE);
  if(root==ALL_ONES)          // all 1's matrix consists of root only, 
    return 1;                 // tree size is 1 no subtree info to check
  
  // we have a non-singleton tree T to traverse 
  // compute number of children
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);
  // read subtree information if available 
  uint64_t *subtree_info = &(z->v[z->n]); // array with subtree_info information
  z->n += nchildren-1;                    // advance z->n to the subtree encoding area
  size_t encode_seen = nchildren-1;       // consume one item x non-last children 
  // visit children 
  size_t cnum = 0;             // current child 
  size_t tree_size = 1;        // account for root node
  for(int i=0;i<4;i++) 
    // invariant: both *pos and z->n point to the beginning of subtree cnum
    if(root & (1<<i)) {
      size_t pc = *pos;  // save current position in m and z 
      size_t nc = z->n;
      size_t child_subtree_size=0, child_encode_size = 0;
      // compute size of subtree encoding
      if(cnum<nchildren-1) {
        child_encode_size = subtree_info[cnum]>>BITSxTSIZE;
        encode_seen += subtree_info[cnum]>>BITSxTSIZE;
      }
      else  // last child 
        child_encode_size = tot_encode_size -encode_seen; // remaining encoding 
      // ------- go down one level ----------- 
      if(size==2*MMsize)  // end of recursion
        *pos += (child_subtree_size = Minimat_node_ratio); // update *pos and child_subtree_size
      else if(child_encode_size==0) {
        k2dfs_visit_fast(size/2,m,pos);  // advance pos to the end of subtree
        child_subtree_size = *pos -pc;   // recover subtree size from advancement in *pos
      }
      else {// recurse on subtree
        child_subtree_size =  k2dfs_check_sizes(size/2,m,pos,z,child_encode_size);
        // check that child_subtree_size matches the advancement in *pos
        if(child_subtree_size != *pos -pc) 
          fprintf(stderr,"Subtree scanned size: %zu, reported size: %zu\n",*pos-pc,child_subtree_size);
      }
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
  (void) pos_save; // avoid warning
  return tree_size;
}
#else
// alternative version in which the subtree size is stored also for the last child
size_t k2dfs_check_sizes(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, 
                                size_t tot_encode_size)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos);        // implies m is non-empty
  size_t pos_save = *pos;     // save starting position of tree
  node_t root = k2read_node(m,*pos); (*pos)++;
  assert(root<ILLEGAL_NODE);
  if(root==ALL_ONES)          // all 1's matrix consists of root only, (this case should not happen)
    return 1;                 // tree size is 1 no subtree info to check
  
  // we have a non-singleton tree T to traverse 
  // compute number of children
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);
  // read subtree information 
  uint64_t *subtree_info = &(z->v[z->n]); // array with subtree_info information
  z->n += nchildren;                      // advance z->n to the subtree encoding area
  size_t encode_seen = nchildren;         // consume one item for each children 
  // visit children 
  size_t cnum = 0;             // current child 
  size_t tree_size = 1;        // account for root node
  for(int i=0;i<4;i++) 
    // invariant: both *pos and z->n point to the beginning of subtree cnum
    if(root & (1<<i)) {
      size_t pc = *pos;  // save current position in m and z 
      size_t nc = z->n;
      size_t child_subtree_size=0;
      // compute size of subtree encoding
      size_t  child_encode_size = subtree_info[cnum]>>BITSxTSIZE;
      encode_seen += child_encode_size;
      // ------- go down one level ----------- 
      if(size==2*MMsize) { // end of recursion
        assert(child_encode_size==0);  // for hieght 2 nodes there canot be a subtree encoding
        *pos += (child_subtree_size = Minimat_node_ratio); // update *pos and child_subtree_size
      } else if(child_encode_size==0) {    // no subtree encoding, scan subtree with dfs
        k2dfs_visit_fast(size/2,m,pos);  // advance pos to the end of subtree
        child_subtree_size = *pos -pc;   // recover subtree size from advancement in *pos
      }
      else {// recurse on subtree
        child_subtree_size =  k2dfs_check_sizes(size/2,m,pos,z,child_encode_size);
        // check that child_subtree_size matches the advancement in *pos
        if(child_subtree_size != *pos -pc) 
          fprintf(stderr,"Subtree scanned size: %zu, reported size: %zu\n",*pos-pc,child_subtree_size);
      }
      // check that stored subtree size matches
      if(child_subtree_size!=(subtree_info[cnum]&TSIZEMASK))
        fprintf(stderr,"Subtree stored size: %zu, reported size: %zu\n",subtree_info[cnum]&TSIZEMASK,child_subtree_size);
      // check stored subtree encoding matches 
      size_t scanned_encoding = z->n-nc;
      if(child_encode_size!=scanned_encoding)
        fprintf(stderr,"Subtree encoding stored size: %zu, scanned size: %zu\n",child_encode_size,scanned_encoding);
      cnum++;
      tree_size += child_subtree_size;
    }  
  assert(cnum==nchildren); // we should have visited all children
  assert(*pos == pos_save + tree_size); // check again tree size
  assert(encode_seen == tot_encode_size); // check that we have seen all subtree encodings
  (void) pos_save; // avoid warning
  (void) tot_encode_size;
  return tree_size;
}
#endif

// similar to k2dfs_check_sizes(), but instead of checking the subtree sizes
// for each node which is destination of a backpointer, we start???? the position
// of each correponding subtree info together with the backpointer 
void k2dfs_compute_backpointer_info(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z)  
{
  #ifdef SIMPLEBACKPOINTERS
  quit("k2dfs_compute_backpointer_info: should not be used for simple backpointers", __LINE__, __FILE__);
  #else
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos);        // implies m is non-empty
  assert(z->n < z->nmax);
  node_t root = k2read_node(m,*pos); 
  assert(root<ILLEGAL_NODE);
  if(*pos > TSIZEMASK) quit("k2dfs_backpointer_info: *pos overflow",__LINE__,__FILE__);
  if(root==ALL_ONES)          // all 1's matrix consists of root only
    return;                   // nothing to do (this case should not happen unless the whole matrix is ALL_ONES)

  // we have a non-singleton tree T with root in *pos and subtring in z->v[z->n] 
  // if pos is destination of a backpointer, the corresponding backpointer
  // is enriched with the value z->n
  pointers_t *ps = m->backp; // backpointer structure
  // skip values in ps->sorted until we reach the first value >= *pos
  while(ps->sidx<ps->size && *pos > ps->nodep[ps->sorted[ps->sidx]]) 
    ps->sidx += 1; // skip all entries smaller than *pos
  while(ps->sidx<ps->size && *pos == ps->nodep[ps->sorted[ps->sidx]]) {
      ps->nodep[ps->sorted[ps->sidx]] |= z->n << BITSxTSIZE; // store the backpointer in z->n
      // uncompress the following line if you want to see the backpointer subtree info on stderr
      //fprintf(stderr,"For node %zu the subtree info is in %zu first size %zu root:%zu, fc %zu\n",*pos, z->n,z->v[z->n]&TSIZEMASK,k2read_node(m,*pos),k2read_node(m,*pos+1));
      ps->sidx += 1;
  }

  // now recurse on the subtrees of T  
  // compute number of children
  (*pos)++;
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);
  // read subtree information 
  uint64_t *subtree_info = &(z->v[z->n]); // array with subtree_info information
  z->n += nchildren;                      // advance z->n to the subtree encoding area
  // visit children 
  size_t cnum = 0;             // current child 
  for(int i=0;i<4;i++) 
    // invariant: both *pos and z->n point to the beginning of subtree cnum
    if(root & (1<<i)) {
      // compute size of subtree encoding
      size_t  child_encode_size = subtree_info[cnum]>>BITSxTSIZE;
      // ------- go down one level ----------- 
      if(size==2*MMsize)              // end of recursion
        *pos += Minimat_node_ratio;    // quickly update *pos 
      else if(child_encode_size==0)     // no subtree encoding, scan subtree with dfs
        k2dfs_visit_fast(size/2,m,pos);    // advance pos to the end of subtree
      else // recurse on subtree
        k2dfs_compute_backpointer_info(size/2,m,pos,z);
      cnum++;
    }  
  assert(cnum==nchildren); // we should have visited all children
  #endif
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

// create and return an interleaved array from the list of entries in a text file
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


// decode to sparse text format a matrix of size 2*MMsize
static void mdecode_to_textfile_base(FILE *outfile, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *a, size_t *pos) 
{
  assert(size==2*MMsize);
  assert(a!=NULL);
  assert(!k2is_empty(a)); 
  assert(!a->main_diag_1 || i==j); // if main_diag_1 is on we must be on a diagonal submatrix
  minimat_t ax[2][2];
  // read root of a
  node_t roota = k2read_node(a,*pos); *pos += 1;
  if(roota==ALL_ONES) {
    if(a->backp!=NULL) quit("Illegal matrix: has backpointers and an ALL_ONES node at last level",__LINE__,__FILE__);
    // output all 1s submatrix  main_diag_1 irrelevant 
    for(size_t ii=i; ii<i+size && ii < msize; ii++)
      for(size_t jj=j; jj<j+size && jj < msize; jj++) {
        int e = fprintf(outfile,"%zu %zu\n",ii,jj);
        if(e<0) quit("mdecode_to_textfile_base: fprintf failed",__LINE__,__FILE__);
      }
  }
  else {
    // split :a taking care also of main_diag
    k2split_minimats(a,pos,roota,ax); // not we cannot pass here i,j
    // fprintf(stderr,"Decoding base submatrix i=%zu j=%zu %x %x %x %x\n",i,j,ax[0][0],ax[0][1],ax[1][0],ax[1][1]);
    for(size_t k=0;k<4;k++) {  
      size_t ii = i + (size/2)*(k/2); size_t jj= j + (size/2)*(k%2);
      minimat_to_text(outfile,msize,ii,jj,size/2, ax[k/2][k%2]); 
    }
  }
}



// recursively decode a k2 submatrix into a list of entries written to a text file
// Parameters:
//   f output file 
//   msize actual file of the matrix
//   i,j submatrix top left corner
//   size k2 submatrix size (has the form 2^k*MMsize)
//   *c input k2mat_t structure
//   *pos position in *c where the submatrix starts
// Note: subtinfo is not used here. if backp!=NULL and a POINTER/ALL_ONES node is found it is followed
// when creating the output if backp==NULL then ALL_ONES nodes are treated as full 1s submatrices
static void mdecode_to_textfile(FILE *outfile, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos) 
{
  assert(!k2is_empty(c)); // never called on an empty matrix
  assert(c->offset==0);   // we are never working on submatrices 
  assert(size%2==0 && size>=2*MMsize);
  assert(i%MMsize==0 && j%MMsize==0);
  assert(i<msize+2*size && j<msize+2*size);
  assert(!c->main_diag_1 || i==j); // if main_diag_1 is on we must be on a diagonal submatrix

  if(size==2*MMsize) { // base case
    mdecode_to_textfile_base(outfile,msize,i,j,size,c,pos);
    return;
  }
  // size>2*MMsize and c contains some data: read c root
  node_t rootc=k2read_node(c,*pos); *pos +=1;
  if(c->backp==NULL && rootc==ALL_ONES) { // all 1s matrix
    // output all 1s submatrix, main_diag_1 irrelevant 
    for(size_t ii=i; ii<i+size && ii < msize; ii++)
      for(size_t jj=j; jj<j+size && jj < msize; jj++) {
        int e = fprintf(outfile,"%zu %zu\n",ii,jj);
        if(e<0) quit("mdecode_to_textfile: fprintf failed",__LINE__,__FILE__);
      }
    return;
  }
  // if pointer node follow it by simply changing position 
  if(c->backp!=NULL && rootc==POINTER) { // recall POINTER=ALL_NODES=0000 
    k2pointer_t destp = k2get_backpointer(c,*pos-1); // -1 because we have already advanced pos
    size_t posp = destp; // move position to the target subtree
    mdecode_to_textfile(outfile,msize,i,j,size,c,&posp);
    return;
  }
  // general case: not a pointer node and there is at least one child
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(size_t k=0;k<4;k++) {  
    size_t ii = i + (size/2)*(k/2); size_t jj= j + (size/2)*(k%2); // submatrix top left corner
    if(rootc & (1<<k)) { // read a submatrix
      if(c->main_diag_1 && (k==0 || k==3)) { // we are on a diagonal submatrix and main_diag_1=true
        k2mat_t tmp = *c;            // update realsize and fullsize 
        tmp.fullsize = tmp.fullsize/2;
        if(k==0) tmp.realsize = (c->realsize > tmp.fullsize) ?  tmp.fullsize :  c->realsize;
        else { 
          assert(k==3);
          assert(c->realsize > tmp.fullsize);
          tmp.realsize = c->realsize - tmp.fullsize;
        }  
        mdecode_to_textfile(outfile,msize,ii,jj,size/2,&tmp,pos);
      }
      else { // off diagonal block of main_diag_1=false
        k2mat_t tmp = *c;
        tmp.main_diag_1 = false;  // do not propagate main_diag_1 outside the diagonal
        mdecode_to_textfile(outfile,msize,ii,jj,size/2,&tmp,pos); // swap ii and jj
      }
    }
    else if(c->main_diag_1 && ii==jj) { // empty subm with main diagonal ones
      for(size_t d=ii; d< ii + (size/2) && d<msize; d++) {
        int e = fprintf(outfile,"%zu %zu\n",d,d);
        if(e<0) quit("mdecode_to_textfile: fprintf failed",__LINE__,__FILE__);
      }
    }
  }
}




// recursively decode a k2 submatrix into a list of entries written to a text file
// only work for plain k2 matrices (no backpointers or main_diagonal)
// Parameters:
//   f output file 
//   msize actual file of the matrix
//   i,j submatrix top left corner
//   size k^2 submatrix size (has the form 2^k*MMsize)
//   *c input k2mat_t structure
//   *pos position in *c where the submatrix starts
// Note: it is currently not used, but it works fine for the results of products
// which are always plain k2 matrices with no backp or main_diag
void mdecode_to_textfile_plain(FILE *outfile, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos) {
  assert(size%2==0 && size>=2*MMsize);
  assert(i%MMsize==0 && j%MMsize==0);
  assert(i<msize+2*size && j<msize+2*size);
  if(c->backp!=NULL) 
    quit("mdecode_to_textfile_plain: cannot handle k2 matrices with backpointers",__LINE__,__FILE__);
  if(c->main_diag_1) 
    quit("mdecode_to_textfile_plain: cannot handle k2 matrices with main_diagonal on",__LINE__,__FILE__); 
  // read c root
  node_t rootc=k2read_node(c,*pos); *pos +=1;
  if(rootc==ALL_ONES) { // all 1s matrix
    for(size_t ii=0; ii<size; ii++)
      for(size_t jj=0; jj<size; jj++)
        if(i+ii<msize && j+jj<msize) { 
          int e = fprintf(outfile,"%zu %zu\n",i+ii,j+jj);
          if(e<0) quit("mdecode_to_textfile_plain: fprintf failed",__LINE__,__FILE__);
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
        mdecode_to_textfile_plain(outfile,msize,ii,jj,size/2,c,pos);
      }
    }
  }
}

// COMPRESSED K2TREE
// later i will do it a more pulished way
// from succint repository
const uint8_t debruijn64_mapping[64] = {
  63,  0, 58,  1, 59, 47, 53,  2,
  60, 39, 48, 27, 54, 33, 42,  3,
  61, 51, 37, 40, 49, 18, 28, 20,
  55, 30, 34, 11, 43, 14, 22,  4,
  62, 57, 46, 52, 38, 26, 32, 41,
  50, 36, 17, 19, 29, 10, 13, 21,
  56, 45, 25, 31, 35, 16,  9, 12,
  44, 24, 15,  8, 23,  7,  6,  5
};

const uint64_t debruijn64 = 0x07EDD5E59A4E28C2ULL;

uint8_t bit_position(uint64_t x){
    return debruijn64_mapping[(x * debruijn64) >> 58];
}

uint8_t _msb(uint64_t x, unsigned long* ret){
  if (!x)
    return false;

  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;

  x ^= x >> 1;
  *ret = bit_position(x);

  return true;
}

uint8_t msb(uint64_t x){
  unsigned long ret = -1U;
  _msb(x, &ret);
  return (uint8_t)ret;
}

uint64_t ceil_log2(uint64_t x) {
  return (x > 1) ? msb(x - 1) + 1 : 0;
}

uint64_t floor_log2(uint64_t x) {
  return (x > 1) ? msb(x) : 0;
}

typedef struct {
  uint64_t n;
  uint64_t* e;
} dsu;

void dsu_init(dsu* u, uint64_t n) {
  u->e = (uint64_t*) malloc(sizeof(int64_t) * n);
  u->n = n;
  for(size_t i = 0; i < n; i++) u->e[i] = i;
}

uint64_t dsu_find_set(dsu* u, uint64_t x) {
  if(u->e[x] == x) return x;
  u->e[x] = dsu_find_set(u, u->e[x]);
  return u->e[x];
}

int8_t dsu_union_set(dsu* u , uint64_t x, uint64_t y) {
  x = dsu_find_set(u, x);
  y = dsu_find_set(u, y);
  if(x == y) return 0;
  if(x > y) {
    int64_t aux = x;
    x = y;
    y = aux;
  }
  u->e[y] = x;
  return 1;
}

void dsu_free(dsu* u) {
  free(u->e);
}

uint64_t get_size_rec(uint8_t *tree, vu64_t *z, uint64_t t_pos, uint64_t t_size,
                      uint64_t z_pos, uint64_t b_tree) {
  uint64_t c_root = __builtin_popcountll(tree[t_pos - 1]);
  if(t_pos == b_tree) {
    if(c_root > 1)
      return z->v[z_pos] & TSIZEMASK;
    else
      return t_size - 1;
  }

  // search for which tree you are
  uint64_t z_m = c_root - 1;
  uint64_t t_m = 0;
  for(uint64_t i = z_pos; i < z_pos + c_root - 1; i++) {
    uint64_t st_size = z->v[i] & TSIZEMASK;
    uint64_t z_skip = z->v[i] >> BITSxTSIZE;

    // found the tree
    if(t_pos + t_m + st_size > b_tree) {
      // move to the first child of this subtree
      return get_size_rec(tree, z, t_pos + t_m + 1, st_size, z_pos + z_m, b_tree);
    }

    // is the brother but not the last
    if(t_pos + t_m + st_size == b_tree && i < z_pos + c_root - 2) {
      return z->v[i + 1] & TSIZEMASK;
    // is the last
    } else if(t_pos + t_m + st_size == b_tree) {
      return t_size - t_m - st_size - 1;
    }
    // moving to next tree
    t_m += st_size;
    z_m += z_skip;
  }
  // is in the last child of a node
  return get_size_rec(tree, z, t_pos + t_m + 1, t_size - t_m - 1, z_pos + z_m, b_tree);
}

static uint64_t get_size(uint8_t *tree, uint64_t t_size, vu64_t *z, uint64_t b_tree) {
  if(b_tree == 0) return t_size;

  uint64_t t_pos = 1;
  uint64_t z_pos = 0;

  return get_size_rec(tree, z, t_pos, t_size, z_pos, b_tree);
}

void k2dfs_write_in_text(size_t size, const k2mat_t *m, size_t *pos, uint8_t* text, size_t *pos_t, uint8_t lvl) {
  assert(size > MMsize);
  assert(size % 2==0);
  assert(*pos < m->pos); // implies m is non-empty
  node_t root = k2read_node(m,*pos); (*pos)++;
  text[*pos_t] = lvl;
  (*pos_t)++;
  text[*pos_t] = (uint8_t) root;
  (*pos_t)++;

  if(root == POINTER) {
    return; // all 1's matrix consists of root only
  }
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      if(size==2*MMsize) { // end of recursion
        minimat_t mm = k2read_minimat(m,pos); // read minimat and advance pos
        text[*pos_t] = lvl + 1;
        (*pos_t)++;
        text[*pos_t] = (uint8_t) mm;
        (*pos_t)++;
      }
      else { // recurse on submatrix
        k2dfs_write_in_text(size/2,m,pos, text, pos_t, lvl + 1); // read submatrix and advance pos
      }
    }
}

// threshold is the minimum amount of nodes to considering erasing a subtree
// block_size is the block size of rank 0000 datastructure
void k2compress(k2mat_t *a, k2mat_t *ca, uint32_t threshold, uint32_t block_size) {
  size_t asize = a->fullsize;
  uint64_t lvs = ceil_log2((uint64_t)asize);
  vu64_t z;
  vu64_init(&z);
  size_t pos = 0;

  // compute all the subtree information 
  k2dfs_sizes(asize, a, &pos, &z, (uint32_t) lvs);

  uint8_t *text = (uint8_t*) malloc(sizeof(uint8_t) * a->pos * 2);
  uint8_t *text2 = (uint8_t*) malloc(sizeof(uint8_t) * a->pos);
  for(size_t i = 0; i < a->pos; i++) {
    text2[i] = (uint8_t) k2read_node(a, i);
  }

  pos = 0;
  size_t pos_t = 0;
  k2dfs_write_in_text(asize, a, &pos, text, &pos_t, 16);
  assert(pos == a->pos);

  int64_t *csa = malloc(sizeof(int64_t) * a->pos * 2);
  int64_t *plcp = malloc(sizeof(int64_t) * a->pos * 2);
  int64_t *lcp = malloc(sizeof(int64_t) * a->pos * 2);

  if(libsais64(text, csa, a->pos * 2, 0, NULL) != 0)
    quit("error creating csa", __LINE__, __FILE__);
  if(libsais64_plcp(text, csa, plcp, a->pos * 2) != 0)
    quit("error creating plcp", __LINE__, __FILE__);
  if(libsais64_lcp(plcp, csa, lcp, a->pos * 2) != 0)
    quit("error creating lcp", __LINE__, __FILE__);

  dsu u;
  dsu_init(&u, a->pos);

  for(size_t i = 1; i < a->pos * 2; i++) {
    uint64_t curr_start_pos = csa[i];
    if(text[curr_start_pos] < 16) {
      continue; // invalid begin
    }

    uint64_t size_subtree = get_size(text2, a->pos, &z, csa[i] / 2);
    // ignoring |trees| <= threshold / 4
    if(size_subtree <= threshold / 4) {
      continue;
    }

    // check that the tree are same length
    if(lcp[i] >= size_subtree * 2) {
      dsu_union_set(&u, curr_start_pos / 2, csa[i - 1] / 2);
    }
  }

  free(csa);
  free(plcp);
  free(lcp);

  uint64_t* prefix_help = (uint64_t*) malloc(sizeof(uint64_t) * a->pos);
  for(size_t i = 0; i < a->pos; i++) prefix_help[i] = 0;

  vu64_t P_h;
  vu64_init(&P_h);

  for(size_t i = 0; i < a->pos; i++) {
    uint64_t repre = dsu_find_set(&u, i);
    if(repre != i) {
      vu64_grow(&P_h, 1);
      P_h.v[P_h.n - 1] = (uint32_t) repre - prefix_help[repre - 1];
      k2add_node(ca, 0);
      size_t next_i = i + get_size(text2, a->pos, &z, i) - 1;
      for(size_t fill = i; fill <= next_i; fill++) {
        prefix_help[fill] = prefix_help[fill - 1] + 1;
      }
      prefix_help[next_i]--;
      i = next_i;
    } else {
      if(i > 0) prefix_help[i] = prefix_help[i - 1];
      k2add_node(ca, text2[i]);
    }
  }
  
  free(prefix_help);
  ca->backp = P_h.n>0 ? pointers_init(&P_h) : NULL;
  vu64_free(&P_h);
  vu64_free(&z);
  free(text);
  dsu_free(&u);
  // no need to compute them, they are not saved
  // rank_init(&(ca->r), block_size, ca);

}

// copy to :a the content of :ca expanding all backpointers 
// work only for full matrices with backpointers (no ALL_ONES nodes)
// used to decompress a matrix just compressed with k2compress
// ignore main_diag_1 
// see k2copy_rec for a version supporting ALL_ONES if backp==NULL
// see k2copy_normalize for a version working for all input matrices 
void k2decompress(size_t size, const k2mat_t *ca, size_t *pos, k2mat_t *a) {
  assert(ca->backp!=NULL);
  assert(ca->r!=NULL);
  assert(!ca->main_diag_1);
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<ca->pos); // implies m is non-empty
  node_t root = k2read_node(ca,*pos); (*pos)++;
  if(root == POINTER) { // is a pointer
    size_t aux = *pos; // remember where to comback

    size_t rp = rank_rank(ca->r, ca, (*pos) - 1);//it is ca[*pos-1]==POINTER
    assert(rp < ca->backp->size);
    *pos = ca->backp->nodep[rp];
    assert(*pos < ca->pos);
    k2decompress(size, ca, pos, a); // read submatrix and advance pos
    
    // moving back to pointer
    *pos = aux;
    return;
  }
  k2add_node(a, root);
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      if(size==2*MMsize) { // end of recursion
        minimat_t mm = k2read_minimat(ca,pos); // read minimat and advance pos
        k2add_minimat(a, mm);
      }
      else { // recurse on submatrix
        k2decompress(size/2, ca, pos, a); // read submatrix and advance pos
      }
    }
}
