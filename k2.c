
// ----------------------------------------------
// global variable storing the size od a mini-matrix 
// ie. the last level of recursion
static const int MMsize = 2; 
// global variable storing how many nodes takes a minimat
// since nodes have 4 bits (not likely to change)
// this amount ot how many nibbles takes a minimat 
static const int Minimat_node_ratio = 1;
#if ((MMsize*MMsize) != 4*Minimat_node_ratio)
#error "check MMsize"
#endif

// node constants (depend on node arity, here 4) 
#define NO_CHILDREN   0x0   // node representing a submatrix of all 0's 
                            // before normalization
#define ALL_ONES      0x0   // node representing a submatrix of all 1's
                            // it is the same as NO_CHILDREN since
                            // after nomralization that code is available 
#define ALL_CHILDREN  0xF   // node which has all children (4), ie a matrix
                            // with all the 4 submatrices non-empty
#define ILLEGAL_NODE  0x10  // illegal node (more than 4 bits)
// minimat constants (depend on size, here are 2x2)
#define MINIMAT0s      0x0   // minimat containing all 0's  
#define MINIMAT1s      0xF   // minimat containing all 1's  

// the type representing a minimat matrix must be an integer
// so that two minimat matrices can be summed with a single
// bitwise OR operation
typedef uint64_t minimat_t; // explict matrix aka minimat (currently 2x2)
// a leaf node must have a number of bits at least equal to the arity
// of the k2-tree (here 4 and it is unlikely to change)
// here we use an uint64_t but only the 4 lower order bits are ever used  
typedef uint64_t node_t;   // non leaf node 

// struct representing a k2-tree: nodes and minimat are stored in a single
// buffer of bytes. offset is used to create a "pointer" to a submatrix
// without copying the buffer during the splitting phase. 
// read_only is currently used only for these
// pointer matrices. By changing b offset can be restricted to the 0/1 values
// so offeset+read_only can be stores in a single byte to save space.
// use read_only also for the input matrices to avoid accidental changes
// or using the const modifier is enough??? 
typedef struct k2mat {
  uint8_t *b;
  size_t pos;   // position where next node is written
  size_t size;  // number of nodes that can be written without rallocating
  size_t offset;// initial pos in b to be skipped (they are not from this matrix)
                // only read only matrices can have a positive offset 
  bool read_only; // if true write and add operations are not allowed
                  // all matrices created by splitting are read only            
} k2mat_t;
// initialize to an empty writable matrix 
#define K2MAT_INITIALIZER {NULL,0,0,0,false}


void quit(const char *msg, int line, char *file);


// ------------------------------------------------------------------- 
// elementary operationson k2mat structures, operating on single items


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
  m->pos = m->size = 0;
}

// nodes are added at the end of a matrix, read back, 
// and sometimes written after a modification 
// (when a node is added we still don't all its children)

// add a node to the current k2_tree
// return the position where the node was stored 
size_t k2add_node(k2mat_t *m, node_t n)
{
  assert(!m->read_only);
  assert(n<ILLEGAL_NODE);
  // make sure there is space
  if(m->pos >= m->size) {
    assert(m->pos ==m->size);
    m->size = 16+2*m->size;          // more than double number of positions
    m->b = realloc(m->b, m->size/2); // each byte stores two positions
    if(m->b==NULL) quit("Unable to enlarge k2-tree",__LINE__,__FILE__);
  }
  assert(m->pos<m->size);
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
  assert(p<m->pos);
  p+=m->offset;
  assert(p<m->size);
  if(p%2==0)
    return m->b[p/2] & 0xF;
  else 
    return  (m->b[p/2] >> 4) & 0xF;
}

// write n at (existing) position p
void k2write_node(k2mat_t *m, size_t p, node_t n)
{
  assert(!m->read_only);
  assert(n<ILLEGAL_NODE);
  assert(p<m->pos);
  if(p%2==0)
    m->b[p/2] = (m->b[p/2]  & 0xF0) | n;
  else 
    m->b[p/2] = (m->b[p/2]  & 0x0F) | (n<<4);
}  


// minimats are only added to the current matrix or read back: never modified

// store a minimatrix in b
// return the starting position where the matrix is stored
size_t k2add_minimat(k2mat_t *b, minimat_t m)
{
  assert(!b->read_only);
  assert(MMsize>7  ||  m< (1UL<<(MMsize*MMsize)));
  // currently a 2x2 matrix takes exactly 4 bit == one node 
  assert(Minimat_node_ratio==1); 
  return k2add_node(b, (node_t) m);
}

// read minimat matrix starting from position p
minimat_t k2read_minimat(const k2mat_t *b, size_t p) {
  assert(Minimat_node_ratio==1);   // otherwise this code must change
  return (minimat_t) k2read_node(b,p);
}

// read all the minimats of a submatrix starting at position *posa
// having already read the root node 
void k2read_minimats(k2mat_t *a,size_t *posa, node_t roota, minimat_t ax[4])
{
  assert(roota!=ALL_ONES);
  for(int i=0;i<4;i++) 
    if(roota & (1<<i)) {
      ax[i] = k2read_minimat(a,*posa);
      *posa += Minimat_node_ratio;
    }
    else ax[i] = MINIMAT0s;
}



// --------------------------------------------------------------
// complex operations on k2mat struct, operating on submatrices



// visit a submatrix from its root (in *pos)
// count the number of nodes and minimatrices visited
// used to split a matrix into 4 submatrices, but also for debugging/info
void k2dfs_visit(int size, const k2mat_t *m, size_t *pos, size_t *nodes, size_t *minimats)
{
  assert(size>=2*MMsize);
  assert(*pos<m->pos);
  node_t root = k2read_node(m,*pos); (*pos)++;
  (*nodes)++;
  if(root==ALL_ONES) return; // all 1's matrix consits of root only
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


// copy the subtree of a starting at *posa to b
// *posa should always point to the next item to be read
// used when summing two matrices and a submatrix is all zeros 
void k2copy_rec(int size, const k2mat_t *a, size_t *posa, k2mat_t *b)
{
  assert(size>MMsize);
  // copy root node
  node_t roota = k2read_node(a,*posa); *posa++;
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

// clone a k2 (sub)matrix of :a starting at position pos
// creating a read only copy :c which is a pointer inside :a
void clone_k2mat(const k2mat_t *a, size_t start, size_t end, k2mat_t *c)
{
  assert(a!=NULL && c!=NULL);
  assert(!k2is_empty(a));
  assert(k2is_empty(c));
  assert(!c->read_only);
  *c = *a; // copy all fields
  c->offset = start;     // actual starting position of c in buffer
  c->pos = end;          // actual ending position of b in buffer 
  c->read_only = true; // b is read only
}


// split the matrix a into 4 submatrices b[0][0], b[0][1], b[1][0], b[1][1]
// the submatrices are "pointers" inside a, so no memory is allocated
void k2split(int size, const k2mat_t *a, k2mat_t b[2][2])
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
  // this code could change if we want to optimize for all 1's matrices
  if(root==ALL_ONES) {
    // if root is all 1's create 4 all 1's submatrices 
    // point to the same root and stop
    for(int i=0;i<2;i++) for(int j=0;j<2;j++)
      k2clone(a, pos-1, pos, &b[i][j]); // ????? check this
    return;
  }
  // root is not all 1's we have the standard structure
  ssize_t next = pos, nodes=0, minimats=0;
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


// NO LONGER USED normalization done in split_and_rec()
// normalize a k2 matrix representation 
// normalize a (sub)matrix stored in :c with the root node
// stored in position :rootpos. 
// :size is the submatrix size  
// must be called only for  size>MMsize
void k2normalize(int size, k2mat_t *c, size_t rootpos)
{
  assert(!c->read_only);
  // MMsize matrices cannot be normalized
  assert(size>MMsize);
  // there should be something else in addition to root
  assert(k2pos(c) > rootpos+1);
  // read current (non normalized) root
  node_t root = k2read_node(c,rootpos);
  // check if matrix is all 0's
  if(root==NO_CHILDREN) {
    assert(k2pos(c)==rootpos+1); // there is only the root
    k2setpos(c,rootpos);         // delete subtree
    return;
  }
  // check if the matrix is all 1's
  if(root==ALL_CHILDREN) {
    if(size==2*MMsize) {  // end of recursion
      bool all_ones = true;
      for(int i=0;i<4 && all_ones; i++) {
        minimat_t m = k2read_minimat(c,rootpos+1+i);
        if(m!=MINIMAT1s) all_ones=false;
      }
    }
    else {
      if(k2pos(c)== rootpos+5) { // all children all 1's
        root = ALL_ONES;
        k2setpos(c,rootpos);    // delete subtree
        k2add_node(c,root);     // write just ALL_ONES root 
      }
    }
  }
  // matrix already normalized: nothing to do 
}
  

// write error message and exit
void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);

  exit(1);
}
