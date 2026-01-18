/* Core routines for handling square binary matrices using k^2 trees 

   type definitions and prototypes of external visible functions

   Matrix dimensions are assumed to be power of 2, of size at least
   MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _K2TYPEDEFS_H
#define _K2TYPEDEFS_H 1

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>

#define SIMPLEBACKPOINTERS  // use simple backpointers, not the full ones
#include "pointers.h" 
#include "rank_0000.h"


// node constants (depend on node arity, here 4 and not likely to change) 
#define NO_CHILDREN   0x0   // node representing a submatrix of all 0's 
                            // before normalization
#define ALL_ONES      0x0   // node representing a submatrix of all 1's
                            // it is the same as NO_CHILDREN since
                            // after normalization that code is available 
#define POINTER       0x0   // node representing a pointer, is is the same as ALL_ONES
                            // since pointers cannot be used with ALL_ONES option
#define ALL_CHILDREN  0xF   // node which has all children (4), ie a matrix
                            // with all the 4 submatrices non-empty
#define ILLEGAL_NODE  0x10  // illegal node (more than 4 bits)

// an internal node must have a number of bits at least equal to the arity
// of the k2-tree (here 4 and it is unlikely to change)
// here we use an uint64_t but only the 4 lower order bits are ever used  
typedef uint64_t node_t;   // non leaf node 



// struct representing a k2-tree: nodes and minimats are stored in a single
// buffer of bytes. offset is used to create a "pointer" to a submatrix
// without copying the buffer during the splitting phase. 
// read_only is currently used only for these pointer matrices. 
// Note: that the size of the k2mat is not stored in the structure:
// it is maintained externally (why is that?)
// Note: public fields are pos, offset, subtinfo, subtinfo_size, backp, r
typedef struct k2mat {
  uint8_t *b;
  size_t lenb;  // number of nodes that can be written in b without reallocating
  // the following two fields are be used to access a submatrix: b[offset,pos-1]
  size_t pos;   // position where next node is written/ position of last node +1
  size_t offset;// initial nodes in b to be skipped (they are not from this matrix)
                // only read only (sub)matrices can have a positive offset
  // the following three fields are used for the subtree size information  
  uint64_t *subtinfoarray; // subtree size encoding if present, or NULL            
  uint64_t *subtinfo;      // start of subtree information for current submatrix 
                           // only read-only (sub)matrices can have subtinfo!=subtinfoarray
  size_t subtinfo_size;    // size (# elements) of subtinfo for current submatrix
  // the following two fields support pointers to subtrees for compressed matrices  
  pointers_t *backp;    // pointers to repeated subtree information
  rank_0000_t *r;       // rank 0000 auxiliary structure
  // flag to denote if the matrix is read only; if true the matrix is a pointer to another matrix
  // and usually defines a submatrix of the other matrix, using offset,pos, and subtinfo 
  // all matrices created by splitting k2split_k2/k2jumpsplit_k2 are read only  
  bool read_only;   // if true write and add operations are not allowed
  bool open_ended;  // pos strictly greater that the end of the matrix: only for read_only matrices      
  // flags to denote if the matrix is transposed or has 1s in the main diagonal 
  bool transpose, main_diag_1;
} k2mat_t;
// initialize to an empty writable matrix 
#define K2MAT_INITIALIZER {NULL,0, 0,0, NULL,NULL,0, NULL,NULL, false, false, false, false}

// maximum allowed size of a k2 matrix
#define MaxMatrixSize (1UL<<40)

// Constants to store size and esizes in a single entry of the subtse array
#define BITSxTSIZE 40
#define TSIZEMASK ( (((uint64_t) 1)<<BITSxTSIZE) -1 )


// float type used for vector elements in matrix-vector multiplication
typedef double vfloat;

// ======== prototypes ===========

// from minimat.c
// init minimatrices: must be called only once with the size of the minimats
// no k2-related function can be called before this one
void minimat_init(int msize);
// return the current minimat size 
uint32_t minimat_size();
// revert the effects of minimat_init and make it possible to call it again
void minimat_reset();

// from k2aux.c
void k2read_subtinfo(k2mat_t *a, const char *infofile);
void k2add_subtinfo_limit(size_t size, k2mat_t *a, size_t limit);
size_t k2treesize(const k2mat_t *m);

// from k2io.c
// save a k2-matrix to file
void msave_to_file(size_t size, size_t asize, const k2mat_t *a, const char *filename);
// load a k2-matrix from file
size_t mload_from_file(size_t *asize, k2mat_t *a, const char *filename);
// write the content of a k2 matrix in a bbm matrix
void mwrite_to_bbm(uint8_t *m, size_t msize, size_t size, const k2mat_t *a);
// read the uncompressed matrix *m of size msize into the k2mat_t structure *a 
size_t mread_from_bbm(uint8_t *m, size_t msize, k2mat_t *a);
// return number of nonzeros in the matrix
size_t mget_nonzeros(size_t asize, const k2mat_t *a);
// write to :file statistics for a k2 matrix :a with an arbitrary :name as identifier
// return number of nonzeros in the matrix
size_t mshow_stats(size_t size, size_t asize, const k2mat_t *a, const char *mname,FILE *file);
// free a k2 matrix
void matrix_free(k2mat_t *m);
// make a read-only copy of a k2 matrix without allocating new memory
void mmake_pointer(const k2mat_t *a, k2mat_t *c);

// from k2ops.c
// check if two k2 compressed matrices :a and :b are equal
// if a==b return -d, where d>0 is the number of levels traversed  
// if a!=b return the level>=0 containing the first difference
int mequals(size_t size, const k2mat_t *a, const k2mat_t *b);
// sum two k2 matrices a and b writing the result to c
// multiplication is done replacing scalar + by logical or 
void msum(size_t size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c);
// multiply two k2 matrices a and b writing the result to c
// multiplication is done replacing scalar */+ by logical and/or 
void mmult(size_t size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c);
// right mutiply a k2 matrix :a by a vector :x writing the result to :y
// :size is the internal size of the k2 matrices (not the size of the vector 
// which can be smaller and in that case :a is padded with zeros)
void mvmult(size_t asize, const k2mat_t *a, size_t size, double *x, double *y, bool clear_y);

// from k2text.c
void mwrite_to_textfile(size_t msize, size_t size, const k2mat_t *a, char *outname);
size_t mread_from_textfile(size_t *msize, k2mat_t *a, char *iname, size_t xsize);
uint64_t k2dfs_sizes(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, int32_t depth2go);
uint64_t k2dfs_sizes_limit(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, size_t limit);
size_t k2dfs_check_sizes(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z, 
                                size_t tot_encode_size);
void k2dfs_compute_backpointer_info(size_t size, const k2mat_t *m, size_t *pos, vu64_t *z);


// compress k2tree
void k2compress(size_t asize, k2mat_t *a, k2mat_t *ca, uint32_t threshold, uint32_t block_size);
void k2decompress(size_t size, const k2mat_t *ca, size_t *pos, k2mat_t *a);

#endif /* _K2TYPDEFS_H */
