/* Core routines for handling square binary matrices using k^2 trees 

   type definitions and prototypes of external visible functions

   Matrix dimensions are assumed to be power of 2, of size at least
   MMSize (minimatrix size), ie, the size of the last level of recursion. 


   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>


// ----------------------------------------------
// global variable storing the size of a mini-matrix 
// ie. the last level of recursion
// possibly this will be a command line parameter
static const int MMsize = 2; 
// global variable storing how many nodes takes a minimat:
// since nodes have 4 bits (not likely to change)
// this amounts to how many nibbles takes a minimat 
// a 2x2 binary matrix takes one nibble, so here the constant is 1
static const int Minimat_node_ratio = 1;
#if ((MMsize*MMsize) != 4*Minimat_node_ratio)
#error "MMsize vs Minimat_node_ratio mismatch"
#endif

// node constants (depend on node arity, here 4 and not likely to change) 
#define NO_CHILDREN   0x0   // node representing a submatrix of all 0's 
                            // before normalization
#define ALL_ONES      0x0   // node representing a submatrix of all 1's
                            // it is the same as NO_CHILDREN since
                            // after normalization that code is available 
#define ALL_CHILDREN  0xF   // node which has all children (4), ie a matrix
                            // with all the 4 submatrices non-empty
#define ILLEGAL_NODE  0x10  // illegal node (more than 4 bits)
// minimat constants (depend on size, here are 2x2)
#define MINIMAT0s      0x0   // minimat containing all 0's  
#define MINIMAT1s      0xF   // minimat containing all 1's  

// the type representing a (temporary) minimat matrix must be an integer
// so that two minimat matrices can be summed with a single
// bitwise OR operation
typedef uint64_t minimat_t; // explict matrix aka minimat (currently 2x2)
// check that minimat_t is large enough to store a minimat
#if ((MMsize*MMsize) > 64 ) // 64=8*sizeof(minimat_t), update if minimat_t changes
#error "MMsize too large for minimat_t"
#endif

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
