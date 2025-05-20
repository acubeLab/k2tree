/* Core routines for handling square binary matrices using k^2 trees 

   This file contains the defintion of the k2mat_t type and the basic
   operations on it, without reference to particular operations
   
   Matrix sizes in these routines (except k2get_k2size) always refer to an 
   internal k2 representation hence are assumed to be power of 2 and at least
   MMSize (minimatrix size), ie, the size of the last level of recursion
   
   The routines in the file support any k2mat size that fits in a 
   size_t variable
   
   This file is included in k2ops.c and should not be compiled separately

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#include "k2.h"
#endif
#include "minimats.c" // includes k2.h bbm.h

#include "pointers.h"
#include "rank_0000.h"

// ------------------------------------------------------------------- 
// elementary operations on k2mat structures, operating on single fields


// return current pos (where the next item will be written 
size_t k2pos(const k2mat_t *m)
{
  return m->pos;
}

// return the size of the tree encoding the (sub)matrix
size_t k2treesize(const k2mat_t *m)
{
  assert(m->pos >= m->offset);
  return m->pos - m->offset;
}


// delete some nodes resetting m->pos to p
// if p==0 this empties the matrix 
void k2setpos(k2mat_t *m, size_t p)
{
  assert(!m->read_only);
  assert(p<m->pos);
  m->pos = p;
}

// check if a matrix is empty
bool k2is_empty(const k2mat_t *m)
{
  return m->pos == 0;
}

// check if the submatrix starting at pos is empty
bool k2submatrix_empty(const k2mat_t *m, size_t pos)
{
  return pos == m->pos;
}

// make m empty, but keep memory  
void k2make_empty(k2mat_t *m)
{
  assert(!m->read_only);
  assert(m->offset==0);
  m->pos = 0;
}

// free mem, *m still reusable if needed 
static void k2_free(k2mat_t *m)
{
  if(m->read_only) // read only matrices are pointers to other matrices
    quit("Illegal operation: freeing a read only k2-matrix",__LINE__,__FILE__); 
  if(m->b!=NULL) free(m->b);
  m->b=NULL;
  m->pos = m->lenb = m->offset = 0;
  if(m->subtinfo!=NULL) {
    free(m->subtinfo); m->subtinfo=NULL; m->subtinfo_size=0;
  }
  if(m->p != NULL)
    pointers_free(m->p);
  if(m->r != NULL)
    rank_free(m->r);
}

// nodes are added at the end of a matrix:
// since when a node is added we still don't know all its children
// once the corresponding subtree is completed nodes are read back 
// and sometimes re-written after a modification 
// (this is why we keep track where the node was written)

// Hence for nodes we have add/read/write operations 


// add a node to the current k2_tree
// return the position where the node was stored 
size_t k2add_node(k2mat_t *m, node_t n)
{
  assert(!m->read_only);
  assert(n<ILLEGAL_NODE);
  assert(m->lenb%2==0);            // #positions must be even
  // make sure there is space
  if(m->pos >= m->lenb) {
    assert(m->pos ==m->lenb);
    m->lenb = 16+2*m->lenb;          // more than double number of positions 
    m->b = realloc(m->b, m->lenb/2); // each byte stores two positions
    if(m->b==NULL) quit("Unable to enlarge k2-tree",__LINE__,__FILE__);
  }
  assert(m->pos<m->lenb);
  // since a node is stored in 4 bits, we store two nodes in a byte
  if(m->pos%2==0)
    m->b[m->pos/2] = (uint8_t) n;    // note: we are writing 0 at m->pos+1, it's ok we are at the end 
  else
    m->b[m->pos/2] = (m->b[m->pos/2] & 0x0F) | (uint8_t) (n<<4);
  // return position where node was stored and advance by 1
  return m->pos++; 
}

// get node at position p 
node_t k2read_node(const k2mat_t *m, size_t p)
{
  p+=m->offset;
  assert(p<m->pos && p<m->lenb);
  if(p%2==0)
    return m->b[p/2] & 0xF;
  else 
    return  (m->b[p/2] >> 4) & 0xF;
}

// write node n at (existing) position p
void k2write_node(k2mat_t *m, size_t p, node_t n)
{
  assert(!m->read_only);
  assert(n<ILLEGAL_NODE);
  assert(p<m->pos && p<m->lenb);
  assert(m->offset==0);   // if m->offset>0 node should be read only 
  if(p%2==0)
    m->b[p/2] = (m->b[p/2]  & 0xF0) | (uint8_t) n;
  else 
    m->b[p/2] = (m->b[p/2]  & 0x0F) | (uint8_t) (n<<4);
}  


// minimats are only added to a matrix or read back: they are never modified
// so we have add/read operations but there is not a k2write_minimat function

// store a minimatrix in b
// return the starting position where the matrix is stored
size_t k2add_minimat(k2mat_t *b, minimat_t m)
{
  assert(!b->read_only);
  assert(MMsize>7  ||  m< (1UL<<(MMsize*MMsize)));
  if (MMsize==2) {
    // a 2x2 matrix takes exactly 4 bit == one node 
    assert(Minimat_node_ratio==1);
    return k2add_node(b, (node_t) m);
  }
  else {
    assert(MMsize==4);
    assert(Minimat_node_ratio==4);
    // save in big endian format 
    size_t first = k2add_node(b, (node_t) ( (m>>12) &0xF));
    for(int i=8;i>=0;i-=4 ) { // 3 more nibbles 
      k2add_node(b, (node_t) ((m>>i)&0xF));
    }
    return first;
  }
}

// read a minimat matrix starting from position *p
// after reading position *p is set at the end of the minimat 
// used by k2split_minimats and k2copy_rec
static minimat_t k2read_minimat(const k2mat_t *b, size_t *p) {
  if(MMsize==2) {
    assert(Minimat_node_ratio==1);
    return (minimat_t) k2read_node(b,(*p)++);
  }
  else {
    assert(MMsize==4);
    assert(Minimat_node_ratio==4);
    minimat_t res = 0;
    for(int i=0;i<Minimat_node_ratio;i++) {
      res <<= 4;
      res |= (minimat_t) k2read_node(b,(*p)++);
    }
    return res;
  }    
}

// split the submatrix :a starting at position *posa into 4 minimats
// and write them to ax[] assuming we have already read the root node :roota 
// we are implicitly assuming we are at the last level of the tree.
// Called by msum_rec, mequals_rec and mmult_base
void k2split_minimats(const k2mat_t *a, size_t *posa, node_t roota, minimat_t ax[2][2])
{
  assert(roota!=ALL_ONES); // currently called with this assumption, could change in the future
  for(int i=0;i<4;i++) 
    if(roota & (1<<i)) ax[i/2][i%2] = k2read_minimat(a,posa); // k2read_minimat advances posa
    else               ax[i/2][i%2] = MINIMAT0s;  // note all 0s minimats are not stored
}



// --------------------------------------------------------------
// complex operations on k2mat struct, operating on submatrices



// visit a size x size non-empty submatrix starting from its root (in *pos)
// the visit is done recursively in depth first order 
// count the number of nodes and minimatrices visited incrementing *nodes and *minimats
// count the number of nonzero entries incrementing *nz
// it is assumed the matrix is not all 0s and that there is a root node so size>MMsize
// used to split a matrix into 4 submatrices, but also for debugging/info
// Note: the number of nonzeros can be incorrect because of overflows 
void k2dfs_visit(size_t size, const k2mat_t *m, size_t *pos, size_t *nodes, size_t *minimats, size_t *nz, size_t *all1)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos); // implies m is non-empty
  node_t root = k2read_node(m,*pos); (*pos)++;
  (*nodes)++;
  if(root == POINTER) { // is a pointer
    uint32_t aux = (uint32_t) *pos; // remember where to comback
    size_t aux_nodes = *nodes; // to not overcount nodes
    size_t aux_minimats = *minimats; // to not overcount minimats

    uint32_t rp = rank_rank(m->r, m, (uint32_t) (*pos) - 1);
    assert(rp < m->p->p_size);
    *pos = m->p->p[rp];
    assert(*pos < m->pos);
    k2dfs_visit(size, m, pos, nodes, minimats, nz, all1); // read submatrix and advance pos
    
    *nodes = aux_nodes; // get back correct stats
    *minimats = aux_minimats; // get back correct stats
    
    // moving back to pointer
    *pos = aux;
    return; // all 1's matrix consists of root only
  }
  // add a macro later
//  if(root==ALL_ONES) {
//    if(size>UINT32_MAX) fprintf(stderr,"Overflow in # of nonzeros: all 1's submatrix of size %zu\n",size);
//    else *nz += size*size; 
//    *all1 +=1;   // found an all 1's matrix 
//    return; // all 1's matrix consists of root only
//  }
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      if(size==2*MMsize) { // end of recursion
        (*minimats)++;
        minimat_t mm = k2read_minimat(m,pos); // read minimat and advance pos
        *nz += (size_t) __builtin_popcountll(mm);
      }
      else { // recurse on submatrix
        k2dfs_visit(size/2,m,pos,nodes,minimats,nz,all1); // read submatrix and advance pos
      }
    }
}

// as above but does not track nodes, minimates and nonzero
// used to split a matrix into 4 submatrices
void k2dfs_visit_fast(size_t size, const k2mat_t *m, size_t *pos)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos); // implies m is non-empty
  node_t root = k2read_node(m,*pos); (*pos)++;
  if(root==ALL_ONES) {
    return; // all 1's matrix consists of root only
  }
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      if(size==2*MMsize) { // end of recursion
        k2read_minimat(m,pos); // read minimat and advance pos
      }
      else { // recurse on submatrix
        k2dfs_visit_fast(size/2,m,pos); // read submatrix and advance pos
      }
    }
}


// copy the subtree of :a starting at *posa to :b
// *posa should always point to the next item to be read
// it is assumed the matrix is not all 0s and that there is a root node so size>MMsize
// used when summing two matrices and a submatrix is all zeros 
void k2copy_rec(size_t size, const k2mat_t *a, size_t *posa, k2mat_t *b)
{
  assert(size>MMsize);
  assert(*posa<a->pos); // implies a is non-empty
  // copy root node
  node_t roota = k2read_node(a,*posa); *posa +=1;
  k2add_node(b,roota);
  if(roota==ALL_ONES) return;  // all 1's matrix consists of root only
  for(int i=0;i<4;i++) {
    if((roota & (1<<i))!=0) {
      if(size==2*MMsize) { // end of recursion
        minimat_t m = k2read_minimat(a,posa);
        k2add_minimat(b,m);
      }
      else { // recurse on submatrix
        k2copy_rec(size/2,a,posa,b);
      }
    }
  }
}

// clone a k2 (sub)matrix of :a starting at position :start and ending at :end-1
// creating a read only copy :c which is a pointer inside :a
// used only by k2split_k2 
static void k2clone(const k2mat_t *a, size_t start, size_t end, k2mat_t *c)
{
  assert(a!=NULL && c!=NULL);
  assert(!k2is_empty(a));
  assert(k2is_empty(c));
  assert(!c->read_only);
  *c = *a; // copy all fields
  c->pos = c->offset + end;     // actual ending position of c in buffer
  c->offset += start;           // actual starting position of c in buffer
  c->subtinfo = NULL;           // if necessary will be initialized later 
  c->subtinfo_size = 0;         // if necessary will be initialized later 
  c->read_only = true;   // c is read only
}

// make c an identical read-only image of matrix a
// the previous content of c is freed and lost 
static void k2make_pointer(const k2mat_t *a, k2mat_t *c)
{
  assert(a!=NULL && c!=NULL);
  k2_free(c);
  *c = *a; // copy all fields (caution, including subtinfo) 
  c->read_only = true;   // c is read only
}

// split the matrix :a into 4 submatrices b[0][0], b[0][1], b[1][0], b[1][1]
// the submatrices are "pointers" inside a, so no memory is allocated
// the submatrices are not minimats, but k2mat_t structs
// since we are not at the last level of the k2 tree
// :a is not all 0's, it could be all 1's (for now)
void k2split_k2(size_t size, const k2mat_t *a, k2mat_t b[2][2])
{
  assert(size>2*MMsize);
  assert(a!=NULL && b!=NULL);
  assert(!k2is_empty(a));
  assert(k2is_empty(&b[0][0]) && k2is_empty(&b[0][1]) &&
         k2is_empty(&b[1][0]) && k2is_empty(&b[1][1]));
  // read root node
  size_t pos = 0;
  node_t root = k2read_node(a,pos); pos++;
  // if root is all 1's create 4 all 1's submatrices and stop
  // this will disappear if we optimize all operations for all 1's submatrices
  if(root==ALL_ONES) {
    // if root is all 1's create 4 all 1's submatrices 
    // pointing to the same root and stop
    for(int i=0;i<2;i++) for(int j=0;j<2;j++)
      k2clone(a, pos-1, pos, &b[i][j]); // ????? check this
    return;
  }
  // root is not all 1's: we have the standard structure
  // and we need to do a real partition of the input matrix 

  // if subtree info is available use it to partition
  // more efficiently and pass the info at the lower levels
  // array of subtree sizes and subtree info
  size_t subt_size[4] = {0,0,0,0};  
  uint64_t *subt_info[4] = {NULL,NULL,NULL,NULL}; 
  size_t subt_info_size[4] = {0,0,0,0};  
  // get number of children
  size_t nchildren = __builtin_popcountll(root);
  assert(nchildren>0 && nchildren<=4);  
  // if we have subtree info fill subt_size[] subt_info[] subt_info_size[] 
  if(a->subtinfo!=NULL) { 
    // start of subtree information
    uint64_t *nextsubtinfo = a->subtinfo + (nchildren-1);
    size_t subt_size_tot = 0, subt_info_size_tot =0;
    for(int i=0;i<nchildren-1;i++) {
      // save subtree size 
      subt_size[i] = a->subtinfo[i] &TSIZEMASK;    // size of subtree i 
      assert(subt_size[i]>0);
      subt_size_tot += subt_size[i];
      // save subtree info size
      subt_info_size[i] = a->subtinfo[i] >>BITSxTSIZE;
      if(subt_info_size[i]>0) {
        subt_info[i] = nextsubtinfo;       // subtree has info, save starting position
        nextsubtinfo += subt_info_size[i]; // start of next subtree info
        subt_info_size_tot += subt_info_size[i];
      }
    }
    // compute size of last subtree
    assert(subt_size_tot +1 < k2treesize(a)); // +1 for the root, there must be a final subtree 
    subt_size[nchildren-1] = k2treesize(a) - subt_size_tot -1; // get size of final subtree
    // consumed information cannot be more than a->subtinfo_size 
    assert(nchildren-1+ subt_info_size_tot <= a->subtinfo_size);
    // remaining information is information for last subtree if any
    subt_info_size[nchildren-1] = a->subtinfo_size - (nchildren-1+ subt_info_size_tot);
    if(subt_info_size[nchildren-1]>0) {
      subt_info[nchildren-1] = nextsubtinfo; // save beginning of last subtree info
      // no need to update nextsubtinfo and subt_info_size_tot
      assert(subt_info[nchildren-1]+subt_info_size[nchildren-1]== a->subtinfo+a->subtinfo_size);
    }
  }
  #ifndef NDEBUG
  size_t nodes=0, minimats=0, nz=0, all1=0;
  #endif
  int child = 0;   // child index, used for subt_size[] subt_info[] 
  size_t next=pos; //now pos==1 since we already read the root
  for(int k=0;k<4;k++) {
    int i=k/2; int j=k%2;
    if(root & (1<<k)) { // k-th child is non empty
      #ifdef NDEBUG
      if(a->subtinfo) next = pos + subt_size[child]; // jump to end of submatrix
      else k2dfs_visit_fast(size/2,a,&next);   // move to end of submatrix
      #else
      k2dfs_visit(size/2,a,&next,&nodes,&minimats,&nz, &all1); // move to end of submatrix
      if(a->subtinfo) assert(next == pos + subt_size[child]);  // check subtree size is as expected
      assert(next==pos+nodes+minimats*Minimat_node_ratio);
      nodes = minimats = 0; // reset number of matrices and nodes
      #endif
      k2clone(a, pos, next, &b[i][j]);        // create pointer to submatrix
      pos = next;                             // advance to next item
      if(a->subtinfo) {
        b[i][j].subtinfo = subt_info[child];    // save subt info if available
        b[i][j].subtinfo_size = subt_info_size[child++];  // save subtinfo size and advance child
      }
    }
  }
  assert(next+a->offset==k2pos(a));
}

// compute the size of the smallest k2mat containing a matrix of size msize
// the size depends on the size of the minimat and grows with powers of 2
// here msize is the actual size of an input matrix that is going to be 
// represented by a k2mat of size the returned value 
size_t k2get_k2size(size_t msize) 
{
  assert(msize>0);
  if(4*Minimat_node_ratio != (MMsize*MMsize)) 
    quit("k2get_size: minimats_init not called",__LINE__,__FILE__);
  size_t s = 2*MMsize;    // size of the smallest legal k2mat
  size_t z = s;           // previous s value to detect unsigned overflow
  while(s < msize) {
    s*=2;
    if(s<z) quit("k2get_k2size: overflow",__LINE__,__FILE__);
    z=s;
  }
  return s;
}

// add the subtree info to a k2 matrix
void k2add_subtinfo(k2mat_t *a, const char *infofile)
{
  assert(infofile!=NULL);
  size_t elsize = sizeof(*(a->subtinfo));  // size of an element in the subtinfo array
  FILE *f = fopen(infofile,"rb");
  if(f==NULL) quit("k2add_subtinfo: cannot open info file", __LINE__,__FILE__);
  if(fseek(f,0,SEEK_END)!=0) quit("k2add_subtinfo: seek error on info file", __LINE__,__FILE__);
  long fsize = ftell(f);
  if(fsize<0) quit("k2add_subtinfo: ftell error on info file", __LINE__,__FILE__);
  if(fsize%elsize!=0) quit("k2add_subtinfo: invalid info file", __LINE__,__FILE__);;
  a->subtinfo_size = fsize/elsize;
  rewind(f);
  a->subtinfo = malloc(fsize);
  if(a->subtinfo==NULL) quit("k2add_subtinfo: cannot allocate memory", __LINE__,__FILE__);
  size_t s = fread(a->subtinfo,elsize,a->subtinfo_size,f);
  if(s!=a->subtinfo_size) quit("k2add_subtinfo: error reading info file", __LINE__,__FILE__);
  fclose(f);
}

// add the backpointers to repeated subtrees in a  k2 matrix
// To be completed
void k2add_backpointers(k2mat_t *a, const char *backfile)
{
  assert(backfile!=NULL);
  // change from here
  size_t elsize = sizeof(*(a->subtinfo));  // size of an element in the subtinfo array
  FILE *f = fopen(backfile,"rb");
  if(f==NULL) quit("k2add_subtinfo: cannot open info file", __LINE__,__FILE__);
  if(fseek(f,0,SEEK_END)!=0) quit("k2add_subtinfo: seek error on info file", __LINE__,__FILE__);
  long fsize = ftell(f);
  if(fsize<0) quit("k2add_subtinfo: ftell error on info file", __LINE__,__FILE__);
  if(fsize%elsize!=0) quit("k2add_subtinfo: invalid info file", __LINE__,__FILE__);;
  a->subtinfo_size = fsize/elsize;
  rewind(f);
  a->subtinfo = malloc(fsize);
  if(a->subtinfo==NULL) quit("k2add_subtinfo: cannot allocate memory", __LINE__,__FILE__);
  size_t s = fread(a->subtinfo,elsize,a->subtinfo_size,f);
  if(s!=a->subtinfo_size) quit("k2add_subtinfo: error reading info file", __LINE__,__FILE__);
  fclose(f);
}


