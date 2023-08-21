/* Core routines for handling square binary matrices using k^2 trees 

   type definitions and prototypes of external visible functions

   Matrix dimensions are assumed to be power of 2, of size at least
   MMSize (minimatrix size), ie, the size of the last level of recursion. 


   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>



// node constants (depend on node arity, here 4 and not likely to change) 
#define NO_CHILDREN   0x0   // node representing a submatrix of all 0's 
                            // before normalization
#define ALL_ONES      0x0   // node representing a submatrix of all 1's
                            // it is the same as NO_CHILDREN since
                            // after normalization that code is available 
#define ALL_CHILDREN  0xF   // node which has all children (4), ie a matrix
                            // with all the 4 submatrices non-empty
#define ILLEGAL_NODE  0x10  // illegal node (more than 4 bits)

// an internal node must have a number of bits at least equal to the arity
// of the k2-tree (here 4 and it is unlikely to change)
// here we use an uint64_t but only the 4 lower order bits are ever used  
typedef uint64_t node_t;   // non leaf node 


// struct representing a k2-tree: nodes and minimat are stored in a single
// buffer of bytes. offset is used to create a "pointer" to a submatrix
// without copying the buffer during the splitting phase. 
// read_only is currently used only for these pointer matrices. 
// By changing b, offset can be restricted to the 0/1 values
// so offeset+read_only could be stored in a single byte to save space.
// ??? use read_only also for the input matrices to avoid accidental changes
// or using the const modifier is enough??? 
typedef struct k2mat {
  uint8_t *b;
  size_t pos;   // position where next node is written
  size_t size;  // number of nodes that can be written without rallocating
  size_t offset;// initial nodes in b to be skipped (they are not from this matrix)
                // only read only matrices can have a positive offset 
  bool read_only; // if true write and add operations are not allowed
                  // all matrices created by splitting are read only            
} k2mat_t;
// initialize to an empty writable matrix 
#define K2MAT_INITIALIZER {NULL,0,0,0,false}


// prototypes

// init minimatrices: must be called once with the size of the minimats
void minimat_init(int msize);
// compute the size of the smallest k2mat containing a matrix of size msize
int k2get_k2size(int msize);
// read the uncompressed matrix *m of size msize into the k2mat_t structure *a 
int mread_uncompressed(uint8_t *m, int msize, k2mat_t *a);
// multiply two k2 matrices a and b writing the result to c
// multiplication is done replacing scalar */+ by logical and/or 
void mmult(int size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c);
// sum two k2 matrices a and b writing the result to c
// multiplication is done replacing scalar + by logical or 
void msum(int size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c);
// check if two k2 compressed matrices :a and :b are equal
// if a==b return -d, where d>0 is the number of levels traversed  
// if a!=b return the level>=0 containing the first difference
int mequals(int size, const k2mat_t *a, const k2mat_t *b);
// get statistics on a matrix
int mstats(int size, const k2mat_t *a, size_t *pos, size_t *nodes, size_t *minimats);

