/* Core routines for handling square binary matrices using dense bit arrays

   type definitions and prototypes of external visible functions
   see b128ops.c for details

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>


typedef __uint128_t uint128_t;


// struct representing a binary matrix with a single bit array
// each matrix row is represented by :colb uint128's
// and the whole matrix by :size*colb int128's
typedef struct b128mat {
  uint128_t *b;   // bit array  
  size_t size;    // size of the matrix
  uint32_t colb;       // # column blocks, ie (size+127)/128
  bool read_only; // if true matrix cannot be overwritten or freed
} b128mat_t;
// initialize to an empty writable matrix 
#define B128MAT_INITIALIZER {NULL,0,0,false}

// maximum allowed size of a bitarray matrix
#define MaxMatrixSize (1UL<<30)

// prototypes

// do nothing, added for compatibility with k2mat
void minimat_init(int);
void minimat_reset(void);
// save a b128-matrix to file
void msave_to_file(const b128mat_t *a, const char *filename);
// load a b128-matrix from file
size_t mload_from_file(b128mat_t *a, const char *filename);
// write the content of a b128 matrix in a bbm matrix
void mwrite_to_bbm(uint8_t *m, const b128mat_t *a);
// read the uncompressed matrix *m of size msize into the b128mat_t structure *a 
size_t mread_from_bbm(uint8_t *m, size_t msize, b128mat_t *a);
// write to :file statistics for a b128 matrix :a with an arbitrary :name as identifier
// return number of nonzeros in the matrix
size_t mshow_stats(const b128mat_t *a, const char *mname,FILE *file);
// check if two b128 compressed matrices :a and :b are equal
// if a==b return -1
// if a!=b return the row index>=0 containing the first difference
int mequals_plain(size_t size, const b128mat_t *a, const b128mat_t *b);
// sum two b128 matrices a and b writing the result to c
// multiplication is done replacing scalar + by logical or 
void msum(const b128mat_t *a, const b128mat_t *b, b128mat_t *c);
// multiply two b128 matrices a and b writing the result to c
// multiplication is done replacing scalar */+ by logical and/or 
void mmult(const b128mat_t *a, const b128mat_t *b, b128mat_t *c);
// free a b128 matrix
void matrix_free(b128mat_t *m);
// make a read-only copy of a b128 matrix without allocating new memory
void mmake_pointer(const b128mat_t *a, b128mat_t *c);
// add identity matrix to a
void madd_identity(b128mat_t *a);
// creates a size x size zero matrix
b128mat_t mat_zero(b128mat_t *b);
// creates a size x size identity matrix
b128mat_t mat_identity(b128mat_t *b);

// from k2text.c
void mwrite_to_textfile(const b128mat_t *a, char *outname);
size_t mread_from_textfile(b128mat_t *a, char *iname, size_t xsize);
