/* Core routines for handling square binary matrices using k^2 trees 

   This file contains the defintion of the k2mat_t type and the basic
   operations on it, without reference to particular operations

   Matrix dimensions are assumed to be power of 2, of size at least
   MMSize (minimatrix size), ie, the size of the last level of recursion. 

   This file can be comiled separately or included as in k2ops.c


   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "minimats.c"



// ------------------------------------------------------------------- 
// elementary operationson on k2mat structures, operating on single items


// return current pos (where the next item will be written 
size_t k2pos(const k2mat_t *m)
{
  return m->pos;
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
  m->pos = 0;
}

// free mem, *m still reusable if needed 
void k2_free(k2mat_t *m)
{
  assert(!m->read_only); // read only matrices are pointers to other matrices
  if(m->b!=NULL) free(m->b);
  m->b=NULL;
  m->pos = m->lenb = 0;
}

// nodes are added at the end of a matrix:
// since when a node is added we still don't all its children
// once the corrsposnding subtree is completed nodes are read back 
// and sometimes re-written after a modification 
// (this is why we keep track where the node was written)

// So for nodes we have add/read/write operations 


// add a node to the current k2_tree
// return the position where the node was stored 
size_t k2add_node(k2mat_t *m, node_t n)
{
  assert(!m->read_only);
  assert(n<ILLEGAL_NODE);
  // make sure there is space
  if(m->pos >= m->lenb) {
    assert(m->pos ==m->lenb);
    m->lenb = 16+2*m->lenb;          // more than double number of positions
    assert(m->lenb%2==0);             // #positions must be even 
    m->b = realloc(m->b, m->lenb/2); // each byte stores two positions
    if(m->b==NULL) quit("Unable to enlarge k2-tree",__LINE__,__FILE__);
  }
  assert(m->pos<m->lenb);
  // since a node is stored in 4 bits, we store two nodes in a byte
  if(m->pos%2==0)
    m->b[m->pos/2] = n;
  else
    m->b[m->pos/2] |= (n<<4);
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
    m->b[p/2] = (m->b[p/2]  & 0xF0) | n;
  else 
    m->b[p/2] = (m->b[p/2]  & 0x0F) | (n<<4);
}  


// minimats are only added to a matrix or read back: they are never modified
// so we have add/read operations but there is not a k2write_minimat function

// store a minimatrix in b
// return the starting position where the matrix is stored
size_t k2add_minimat(k2mat_t *b, minimat_t m)
{
  assert(!b->read_only);
  assert(MMsize>7  ||  m< (1UL<<(MMsize*MMsize)));
  // currently a 2x2 matrix takes exactly 4 bit == one node 
  assert(Minimat_node_ratio==1); // otherwise this code must change
  return k2add_node(b, (node_t) m);
}

// read minimat matrix starting from position p
// used by k2read_minimats and k2copy_rec
static minimat_t k2read_minimat(const k2mat_t *b, size_t p) {
  assert(Minimat_node_ratio==1);   // otherwise this code must change
  return (minimat_t) k2read_node(b,p);
}

// split the submatrix :a starting at position *posa into 4 minimats
// and write them to ax[] assuming we have already read the root node :roota 
// we are implicitly assuming we are at the last level of the tree
// called by msum_rec, mequals_rec and mmult_base
void k2split_minimats(const k2mat_t *a, size_t *posa, node_t roota, minimat_t ax[2][2])
{
  assert(roota!=ALL_ONES); // currently called with this assumption, could change in the future
  for(int i=0;i<4;i++) 
    if(roota & (1<<i)) {
      ax[i/2][i%2] = k2read_minimat(a,*posa);
      *posa += Minimat_node_ratio;
    }
    else ax[i/2][i%2] = MINIMAT0s;
}



// --------------------------------------------------------------
// complex operations on k2mat struct, operating on submatrices



// visit a size x size non-empty submatrix starting from its root (in *pos)
// the visit is done recursively in depth first order 
// count the number of nodes and minimatrices visited incrementing *nodes and *minimats
// it is assumed the matrix is not all 0s and that there is a root node so size>MMsize
// used to split a matrix into 4 submatrices, but also for debugging/info
void k2dfs_visit(int size, const k2mat_t *m, size_t *pos, size_t *nodes, size_t *minimats)
{
  assert(size>MMsize);
  assert(*pos<m->pos); // implies m is non-empty
  node_t root = k2read_node(m,*pos); (*pos)++;
  (*nodes)++;
  if(root==ALL_ONES) return; // all 1's matrix consists of root only
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      if(size==2*MMsize) { // end of recursion
        (*minimats)++;
        *pos += Minimat_node_ratio;
      }
      else { // recurse on submatrix
        k2dfs_visit(size/2,m,pos,nodes,minimats);
      }
    }
}


// copy the subtree of :a starting at *posa to :b
// *posa should always point to the next item to be read
// it is assumed the matrix is not all 0s and that there is a root node so size>MMsize
// used when summing two matrices and a submatrix is all zeros 
void k2copy_rec(int size, const k2mat_t *a, size_t *posa, k2mat_t *b)
{
  assert(size>MMsize);
  assert(*posa<a->pos); // implies a is non-empty
  // copy root node
  node_t roota = k2read_node(a,*posa); *posa +=1;
  k2add_node(b,roota);
  if(roota==ALL_ONES) return;  // all 1's matrix consits of root only
  for(int i=0;i<4;i++) {
    if((roota & (1<<i))!=0) {
      if(size==2*MMsize) { // end of recursion
        minimat_t m = k2read_minimat(a,*posa);
        *posa += Minimat_node_ratio;
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
  c->offset = start;     // actual starting position of c in buffer
  c->pos = end;          // actual ending position of c in buffer 
  c->read_only = true;   // c is read only
}


// split the matrix :a into 4 submatrices b[0][0], b[0][1], b[1][0], b[1][1]
// the submatrices are "pointers" inside a, so no memory is allocated
// the submatrices are not minimats, but k2mat_t structs
// so we are not at the last level of the tree
// :a is not all 0's, it could be all 1's (for now)
void k2split_k2(int size, const k2mat_t *a, k2mat_t b[2][2])
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
  // this will disappear if we optimize all operations for all 1's matrices
  if(root==ALL_ONES) {
    // if root is all 1's create 4 all 1's submatrices 
    // pointing to the same root and stop
    for(int i=0;i<2;i++) for(int j=0;j<2;j++)
      k2clone(a, pos-1, pos, &b[i][j]); // ????? check this
    return;
  }
  // root is not all 1's: we have the standard structure
  size_t next = pos, nodes=0, minimats=0;
  for(int k=0;k<4;k++) {
    int i=k/2; int j=k%2;
    if(root & (1<<k)) { // k-th child is non empty
      k2dfs_visit(size/2,a,&next,&nodes,&minimats); // get size of submatrix
      assert(next==pos+nodes+minimats*Minimat_node_ratio);
      k2clone(a, pos, next, &b[i][j]);         // create pointer to submatrix
      pos = next; // advance to next item
      nodes = minimats = 0; // reset number of matrices and nodes
    }
  }
  assert(next==k2pos(a));
}
