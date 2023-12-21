/* Routines for arithmetic/logical operations on binary matrices 
   represented as k^2 trees 

   This file contains the defintion of the operations that make use of
   the basic operations defined in k2aux.c

   Matrix dimensions are assumed to be power of 2, of size at least
   2*MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "k2aux.c"
#include "bbm.h"



// write the content of the :size x :size k2 matrix :a to the bbm matrix :m 
// of size msize*msize. It is assumed m was already correctly initialized and allocated
void mwrite_to_bbm(uint8_t *m, int msize, int size, const k2mat_t *a);

// compress the matrix *m of size msize into the k2mat_t structure *a 
// m should be an array of size msize*msize 
// the old content of :a is lost
// return the size of the k2 matrix (which has the form 2**k*MMsize)
int mread_from_bbm(uint8_t *m, int msize, k2mat_t *a);


// write to :file statistics for a k2 matrix :a with an arbitrary :name as identifier
void mshow_stats(size_t size, int asize, const k2mat_t *a, const char *mname,FILE *file);



// main entry point for matrix equality
// check if size x size k2 compressed matrices :a and :b are equal
// if a==b return -d, where d>0 is the number of levels traversed  
// if a!=b return the level>=0 containing the first difference
// (first in the sense of the first level encountered in dfs order)
// Note that if a==b we return the number of visited levels (tree depth), 
// while if a!=b we return the level of the first difference counting from 0 (root)
// these two values differ by one (these can be seen in the returned value)
// :a and :b must be of size at least 2*MMsize but their content can be
// arbitrary: all 0's, all 1's, or generic
// note: here all 0's matrices are considered of depth 1 even if they are empty
int mequals(int size, const k2mat_t *a, const k2mat_t *b);



// main entry point for matrix addition
// sum size x size k2 compressed matrices :a and :b storing
// the result in :c, the old content of :c is discarded
// :a and :b must be of (same) size at least 2*MMsize but their content can be 
// arbitrary: all 0's, all 1's, or generic
// at exit:
//    if the result is a zero matrix c is left empty
//    if the result is an all one's matrix c contains a single ALL_ONES node
//    otherwise c is a node + the recursive description of its subtree  
void msum(int size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c);


// main entry point for matrix multiplication. 
// multiply size x size k2 compressed matrices :a and :b storing
// the result in :c, old content of :c is discarded
// :a and :b must be of size at least 2*MMsize but their content can be 
// arbitrary: all 0's, all 1's, or generic
// at exit:
//    if the result is a zero matrix c is left empty
//    if the result is an all one's matrix c contains a single ALL_ONES node
//    otherwise c is a node + the recursive description of its subtree  
void mmult(int size, const k2mat_t *a, const k2mat_t *b, k2mat_t *c);


// save the matrix :a to file :filename
// the format is
// 1) the actual size of the matrix (an int)
// 2) the size of the last level of recursion, aka the size of a minimat,  (an int)
// 3) the size of the k2 matrix (an int) 
// 4) the total number of positions (a size_t)
// 5) the array of positions (an uint8_t array, each uint8 stores 2 positions)
void msave_to_file(int size, int asize, const k2mat_t *a, const char *filename);

// load a k2 matrix stored in file :filename into the k2mat_t structure :a
// return the actual size of the matrix, see msave_to_file() for the file format
// :a old content is discarded
// must be called after minimat_init() since it checks that the correct MMsize is used
int mload_from_file(int *asize, k2mat_t *a, const char *filename);



// non-k2 names to two useful functions

// free mem, *m still reusable if needed 
void matrix_free(k2mat_t *m);

// make c an identical read-only image of matrix a
// the previous content of c is freed and lost 
void mmake_pointer(const k2mat_t *a, k2mat_t *c);


