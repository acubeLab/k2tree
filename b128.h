/* Core routines for handling square binary matrices using k^2 trees 

   type definitions and prototypes of external visible functions

   Matrix dimensions are assumed to be power of 2, of size at least
   MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Technical note: matrix sizes are int, therefore limited to 2^31, but values related to
   the overall number of elements is always stored into a size_t variable (usually 64 bits))  

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>


typedef __uint128_t uint128_t;


// struct representing a b128-tree: nodes and minimat are stored in a single
// buffer of bytes. offset is used to create a "pointer" to a submatrix
// without copying the buffer during the splitting phase. 
// read_only is currently used only for these pointer matrices. 
// By changing b, offset can be restricted to the 0/1 values
// so offeset+read_only could be stored in a single byte to save space.
// ??? use read_only also for the input matrices to avoid accidental changes
// or using the const modifier is enough??? 
typedef struct b128mat {
  uint128_t *b;   // bit array  
  size_t size;    // size of the matrix
  int colb;       // # column blocks, ie (size+127)/128
  bool read_only; // if true write and add operations are not allowed
} b128mat_t;
// initialize to an empty writable matrix 
#define B128MAT_INITIALIZER {NULL,0,0,false}


// prototypes

// save a b128-matrix to file
void msave_to_file(int size, int asize, const b128mat_t *a, const char *filename);
// load a b128-matrix from file
int mload_from_file(int *asize, b128mat_t *a, const char *filename);
// write the content of a b128 matrix in a bbm matrix
void mwrite_to_bbm(uint8_t *m, int msize, int size, const b128mat_t *a);
// read the uncompressed matrix *m of size msize into the b128mat_t structure *a 
int mread_from_bbm(uint8_t *m, int msize, b128mat_t *a);
// write to :file statistics for a b128 matrix :a with an arbitrary :name as identifier
void mshow_stats(size_t size, int asize, const b128mat_t *a, const char *mname,FILE *file);
// multiply two b128 matrices a and b writing the result to c
// multiplication is done replacing scalar */+ by logical and/or 
void mmult(int size, const b128mat_t *a, const b128mat_t *b, b128mat_t *c);
// sum two b128 matrices a and b writing the result to c
// multiplication is done replacing scalar + by logical or 
void msum(int size, const b128mat_t *a, const b128mat_t *b, b128mat_t *c);
// check if two b128 compressed matrices :a and :b are equal
// if a==b return -d, where d>0 is the number of levels traversed  
// if a!=b return the position>=0 containing the first difference
int mequals(int size, const b128mat_t *a, const b128mat_t *b);
// get statistics on a matrix
// int mstats(int size, const b128mat_t *a, size_t *pos, size_t *nodes, size_t *minimats);
// free a b128 matrix
void b128_free(b128mat_t *m);
// make a read-only copy of a b128 matrix without allocating new memory
void b128make_pointer(const b128mat_t *a, b128mat_t *c);

