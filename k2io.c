/* Routines for input/output (and related) operations on binary matrices 
   represented as k^2 trees 

   This file contains the definitions of complex operations that make 
   use of the basic operations defined in k2aux.c

   Matrix dimensions are assumed to be power of 2, of size at least
   2*MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#pragma GCC target ("sse4.2")  // ensure sse4.2 compiler switch it is used 
#include "minimats.c"    // includes minimats.c k2.h bbm.h

static void mdecode_bbm(uint8_t *m, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos);
static void mencode_bbm(uint8_t *m, size_t msize, size_t i, size_t j, size_t size, k2mat_t *c);



// BBM support will be removed in future versions
// write the content of the :size x :size k2 matrix :a to the bbm matrix :m 
// of size msize*msize. It is assumed m was already correctly initialized and allocated
void mwrite_to_bbm(uint8_t *m, size_t msize, size_t size, const k2mat_t *a)
{
  assert(size>=msize);
  if(k2is_empty(a)) {  // an empty k2 matrix is all 0s
    byte_to_bbm(m,msize,0,0,msize,0); // fill m with 0s
    return;
  }
  byte_to_bbm(m,msize,0,0,msize,2); // fill m with illegal value 2
  size_t pos = 0;
  mdecode_bbm(m,msize,0,0,size,a,&pos);
  assert(pos==k2pos(a)); // check we read all the k2mat_t structure
}


// compress the matrix *m of size :msize into the k2mat_t structure *a 
// m should be an array of size msize*msize 
// the old content of :a is lost
// return the size of the k2 matrix (which has the form 2**k*MMsize)
size_t mread_from_bbm(uint8_t *m, size_t msize, k2mat_t *a)
{
  assert(a!=NULL && m!=NULL);
  k2_free(a); // free previous content of a
  assert(msize>1);
  size_t asize = k2get_k2size(msize);
  assert(asize>=2*MMsize);
  // read matrix m into a
  mencode_bbm(m,msize,0,0,asize,a);
  a->fullsize = asize;
  a->realsize = msize;
  return asize;
}

// return statistics on matrix a
// write number of used pos,nodes, minimats and nonzeros in the variables passed by reference
// and return the number of levels
static int mstats(size_t asize, const k2mat_t *a, size_t *pos, size_t *nodes, size_t *minimats, size_t *nz, size_t *all1)
{
  *pos=*nodes=*minimats=*nz=*all1=0;
  if(!k2is_empty(a)) k2dfs_visit(asize,a,pos,nodes,minimats,nz,all1);
   return k2tree_levels(asize,a); // numer of level in tree
}

// return number of nonzero elements in matrix :a
size_t mget_nonzeros(size_t asize, const k2mat_t *a) {
  size_t pos, nodes, minimats, nz, all1;
  mstats(asize,a,&pos,&nodes,&minimats,&nz,&all1);
  return nz;
}

// write to :file statistics for a k2 matrix :a with an arbitrary :name as identifier
// :size is the actual; matrix size (not power of 2), :asize is the internal size
// return number of nonzeros in the matrix
size_t mshow_stats(const k2mat_t *a, const char *mname,FILE *file) {
  size_t pos, nodes, minimats, nz, all1;
  size_t asize=a->fullsize;
  size_t size=a->realsize;
  fprintf(file,"%s:\n matrix size: %zu, leaf size: %d, k2_internal_size: %zu\n",mname,size,MMsize,asize);  
  if(a->main_diag_1) 
    fprintf(file,"matrix has all 1s on the main diagonal: the reported number of nonzeros is a lower bound\n");
  int levels = mstats(asize,a,&pos,&nodes,&minimats,&nz,&all1);
  if(a->backp == NULL)
    assert(pos==nodes+minimats*Minimat_node_ratio); // check that the number of positions is correct
  fprintf(file," #Nonzeros: %zu, Nonzero x row: %lf\n", nz, (double) nz/(double)size);
  fprintf(file," Levels: %d, Nodes: %zu, Minimats: %zu, 1's submats: %zu\n",
          levels,nodes,minimats, all1);
  fprintf(file," Subtree info size (bytes): %zu,", a->subtinfo_size);
  size_t bits_sub = sizeof(*(a->subtinfo)) * a->subtinfo_size;
  fprintf(file, " Subtree info size (bits): %zu\n", bits_sub);
  fprintf(file," #Subtree pointers: %zu\n", a->backp ? a->backp->size : 0);
  size_t bits_p = pointers_size_in_bits(a->backp);
  size_t bits_r = rank_size_in_bits(a->r);
  fprintf(file, " Pointers size (bits): %zu, Rank DS (bits, not stored): %zu\n", bits_p, bits_r);
  // each pos takes 4 bits, so tree size in bytes is (pos+1)/2         
  fprintf(file," Tree size: %zu bytes, %zu bits, Bits x nonzero: %lf\n",
          (pos+1)/2 , 4*pos, 4.0*(double)(pos)/(double) nz);
  fprintf(file, " Total space (bits): %zu, Bits x nonzero: %lf\n", 
      pos * 4 + bits_r + bits_p, (double) (pos * 4 + bits_r + bits_p) / (double) nz);
  return nz;
}

// save the matrix :a to file :filename
// the format is
// 1) the actual size of the matrix (a size_t)
// 2) the size of the last level of recursion, aka the size of a minimat, (an int)
// 3) the internal size of the k2 matrix (a size_t, we could skip this) 
// 4) the total number of positions (a size_t)
// 5) the array of positions (an uint8_t array, each uint8 stores 2 positions)
void msave_to_file(const k2mat_t *a, const char *filename)
{
  assert(a!=NULL);
  size_t size = a->realsize;
  size_t asize = a->fullsize;
  FILE *f = fopen(filename,"w");
  if(f==NULL) quit("msave_to_file: cannot open file", __LINE__,__FILE__);
  size_t e = fwrite(&size,sizeof(size),1,f);
  if(e!=1) quit("msave_to_file: cannot write size",__LINE__,__FILE__);
  e = fwrite(&MMsize, sizeof(MMsize),1,f);
  if(e!=1) quit("msave_to_file: cannot write MMsize",__LINE__,__FILE__);
  e = fwrite(&asize, sizeof(asize),1,f);
  if(e!=1) quit("msave_to_file: cannot write asize",__LINE__,__FILE__);
  e = fwrite(&a->pos,sizeof(size_t),1,f);
  if(e!=1) quit("msave_to_file: cannot write number of positions",__LINE__,__FILE__);
  if(a->pos>0) {
    size_t bytes = (a->pos+1)/2;
    e = fwrite(a->b,sizeof(uint8_t),bytes,f);
    if(e!=bytes) quit("msave_to_file: cannot write positions",__LINE__,__FILE__);
  }
  fclose(f);
}


// load a k2 matrix stored in file :filename into the k2mat_t structure :a
// return the actual size of the matrix and store to *asize the internal (power of 2)
// size of the k2 matrix, see msave_to_file() for details of the file format
// :a old content is discarded
// call minimat_init() if necessary, check that the correct MMsize is used
size_t mload_from_file(k2mat_t *a, const char *filename)
{
  assert(a!=NULL);
  k2_free(a);
  FILE *f = fopen(filename,"r");
  if(f==NULL) quit("mload_from_file: cannot open file", __LINE__,__FILE__);
  size_t size, asize;
  int mmsize;
  size_t e = fread(&size,sizeof(size),1,f);
  if(e!=1) quit("mload_from_file: cannot read matrix size",__LINE__,__FILE__);
  if(size<=1) quit("mload_from_file: matrix size smaller than 2, wrong format?",__LINE__,__FILE__);
  e = fread(&mmsize, sizeof(int),1,f);
  if(e!=1) quit("mload_from_file: cannot read minimatrix size",__LINE__,__FILE__);
  if(MMsize==0)
    minimat_init(mmsize); // initialize minimat library if not already done
  else 
    if(mmsize!=MMsize) quit("mload_from_file: minimatrix size mismatch",__LINE__,__FILE__);
  e = fread(&asize, sizeof(asize),1,f);
  if(e!=1) quit("mload_from_file: cannot read k2 matrix size",__LINE__,__FILE__);
  if(asize<2*MMsize)
    quit("mload_from_file: k2 matrix size incompatible with minimatrix size, wrong format?",
                         __LINE__,__FILE__);
  if(asize!=k2get_k2size(size)) quit("mload_from_file: wrong k2 matrix size",__LINE__,__FILE__ ); 
  a->fullsize = asize; a->realsize = size;                      
  e = fread(&a->pos,sizeof(size_t),1,f);
  if(e!=1) quit("mload_from_file: cannot read number of positions",__LINE__,__FILE__);
  if(a->pos>0) {
    size_t bytes = (a->pos+1)/2;
    a->lenb = 2*bytes;  // equivalent to:   a->lenb = a->pos%2 ? a->pos+1: a->pos;
    a->b = malloc(bytes);
    if(a->b==NULL) quit("mload_from_file: cannot allocate memory",__LINE__,__FILE__);
    e = fread(a->b,sizeof(uint8_t),bytes,f);
    if(e!=bytes) quit("mload_from_file: cannot read positions",__LINE__,__FILE__);
  }
  fclose(f);
  return size;
}


// load a k2 matrix, stored in file :fname, into the k2mat_t structure :a
// possibly read also the subtree info from file :subtname and the backpointers from file :backpname
// if backpname is present, rank_block_size must be provided as it is used to initialize the rank structure
// the old content of :a is discarded
// return the actual size of the matrix and store to *asize the internal (power of 2)
// size of the k2 matrix, see msave_to_file() for details of the file format
size_t mload_extended(k2mat_t *a, char *fname, char *subtname, const char *backpname, uint32_t rank_block_size)
{
  size_t size = mload_from_file(a,fname);
  // try to load subtinfo if present
  if(subtname!=NULL)  k2read_subtinfo(a,subtname);
  if(backpname!=NULL) {
    a->backp = pointers_load_from_file(backpname);
    assert(rank_block_size>0 && rank_block_size%4==0);  
    rank_init(&(a->r),rank_block_size,a);
  }
  return size;
}

// give non-k2 names to two useful functions
// for compatibility with bitarray representation

// free mem, *m still reusable if needed 
void matrix_free(k2mat_t *m) {
  k2_free(m);
}

// make c an identical read-only image of matrix a
// the previous content of c is freed and lost 
void mmake_pointer(const k2mat_t *a, k2mat_t *c)
{
  k2make_pointer(a,c);
}



// ----------- static auxiliary functions ------------

// recursively encode a binary submatrix m[i,i+size)[j,j+size) given in one-byte 
// per entry format into a k2mat_t structure
// Parameters:
//   m matrix to encode of size msize*msize
//   i,j submatrix top left corner
//   size submatrix size (has the form 2^k*MMsize)
//   *c output k2mat_t structure to be filled in dfs order 
static void mencode_bbm(uint8_t *m, size_t msize, size_t i, size_t j, size_t size, k2mat_t *c) {
  assert(size%2==0 && size>=2*MMsize);
  assert(i%MMsize==0 && j%MMsize==0);
  assert(i<msize+2*size && j<msize+2*size);
  // if we are outside m it's an all 0 submatrix and there is nothing to do
  if(i>=msize || j>=msize) return;
  // start building c
  size_t rootpos = k2add_node(c,ALL_ONES);  // write ALL_ONES as root placeholder 
  node_t rootc=NO_CHILDREN;                 // actual root node to be computed
  bool all_ones=true;  // true if all c submatrices cx[i][j] are all 1's
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(int k=0;k<4;k++) {  
    size_t ii = i + (size/2)*(k/2); size_t jj= j + (size/2)*(k%2);
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
      mencode_bbm(m,msize,ii,jj,size/2,c);
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
  else if(all_ones && Use_all_ones_node) {   // all 1s matrix is represented by the ALL_ONES root only
    assert(rootc==ALL_CHILDREN);
    k2setpos(c,rootpos+1);       // discard children
    assert(k2read_node(c,rootpos)==ALL_ONES); // ALL_ONES was the default root 
  }
  else k2write_node(c,rootpos,rootc); // just fix root 
}

// recursively decode a k2 submatrix into a binary submatrix
// m[i,i+size)[j,j+size) in one-byte per entry format 
// Parameters:
//   m output matrix of (overall) size msize*msize
//   i,j submatrix top left corner
//   size submatrix size (has the form 2^k*MMsize)
//   *c input k2mat_t structure
//   *pos position in *c where the submatrix starts
static void mdecode_bbm(uint8_t *m, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos) {
  assert(size%2==0 && size>=2*MMsize);
  assert(i%MMsize==0 && j%MMsize==0);
  assert(i<msize+2*size && j<msize+2*size);
  // read c root
  node_t rootc=k2read_node(c,*pos); *pos +=1;
  if(rootc==ALL_ONES) { // all 1s matrix
    byte_to_bbm(m,msize,i,j,size,1);
    return;
  }
  // here we are assuming that the submatrices are in the order 00,01,10,11
  for(size_t k=0;k<4;k++) {  
    size_t ii = i + (size/2)*(k/2); size_t jj= j + (size/2)*(k%2);
    if(rootc & (1<<k)) { // read a submatrix
      if(size==2*MMsize) { // read a minimatrix
        minimat_t cx = k2read_minimat(c,pos);
        assert(cx!=MINIMAT0s); // should not happen
        minimat_to_bbm(m,msize,ii,jj,size/2,cx);
      }
      else { // decode submatrix
        mdecode_bbm(m,msize,ii,jj,size/2,c,pos);
      }
    }
    else { // the k-th chilldren is 0: write a 0 submatrix
      // if m initialized with 0s we can skip the following line and save some writes
      byte_to_bbm(m,msize,ii,jj,size/2,0);
    }
  }
}

