/* Routines for converting binary matrices in text form, ie 
   a list of directed arcs represented as a pair of vertices,
   to the compressed k^2 tree representation 

   Matrix dimensions are assumed to be power of 2, of size at least
   2*MMSize (minimatrix size), ie, the size of the last level of recursion. 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <inttypes.h>


// prototypes of static functions
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize);
static size_t mread_from_ia(uint64_t ia[], size_t n, size_t msize, k2mat_t *a);
static void mencode_ia(uint64_t *ia, size_t n, uint64_t imin, uint64_t imax, size_t size, k2mat_t *c);
static void mdecode_to_textfile(FILE *outfile, size_t msize, size_t i, size_t j, size_t size, const k2mat_t *c, size_t *pos);
static size_t binsearch(uint64_t *ia, size_t n, uint64_t x);

// read a matrix from the text file :iname (one arc per line)
// and store it in k2 format
// the compressed matrix is stored to :a and it size to :msize
// return the size of the k2 matrix (which has the form 2**k*MMsize)
size_t mread_from_textfile(size_t *msize, k2mat_t *a, char *iname)
{
  assert(iname!=NULL && a!=NULL);
  FILE *f = fopen(iname,"rt");
  if(f==NULL) quit("mread_from_file: cannot open input file",__LINE__,__FILE__);
  // generate interleaved array from input file
  size_t n; // number of arcs
  uint64_t *ia = create_ia(f,&n,msize);
  fclose(f);
  // compress the matrix represented by the ia[] array into a k2mat_t structure
  size_t asize = mread_from_ia(ia,n,*msize,a);
  free(ia);
  return asize;
}



// write the content of the :size x :size k2 matrix :a to the bbm matrix :m 
// of size msize*msize. It is assumed m was already correctly initialized and allocated
void text2bintextfile(size_t msize, size_t size, const k2mat_t *a, char *outname)
{
  assert(outname!=NULL && a!=NULL);
  assert(size>=msize);
  FILE *f = fopen(outname,"wt");
  if(f==NULL) quit("mwrite_to_file: cannot open output file",__LINE__,__FILE__);

  if(k2is_empty(a)) {  // an empty k2 matrix has no arcs
    fclose(f);
    return;
  }
  size_t pos = 0;
  mdecode_to_textfile(f,msize,0,0,size,a,&pos);
  fclose(f);
  assert(pos==k2pos(a)); // check we read all the k2mat_t structure
}

// ----------- static auxiliary functions ------------

// compress the matrix of size msize represented by the interleaved
// array ia[0..n-1] into the k2mat_t structure *a 
// ia[] should be an interleaved array of length n
// the old content of :a is lost
// return the size of the k2 matrix (which has the form 2**k*MMsize)
static size_t mread_from_ia(uint64_t ia[], size_t n, size_t msize, k2mat_t *a)
{
  assert(ia!=NULL && a!=NULL);
  assert(n>0);                // we cannot represent an empty matrix
  assert(msize>1);
  assert( n <= msize*msize); // entries must be less that msize**2
  k2_free(a);                // free previous content of a 
  if(msize>(1<<30)) quit("mread_from_ia: matrix too large",__LINE__,__FILE__);
  size_t asize = k2get_k2size(msize);
  assert(asize>=2*MMsize);
  // encode ia[] into a k2mat_t structure
  mencode_ia(ia,n,0,asize*asize,asize,a);
  return asize;
}


// in a sorted uint64_t array ia[0,n-1] containing distinct values find 
// the first entry >= x using binary search 
// return n if no such entry exists
static size_t binsearch(uint64_t *ia, size_t n, uint64_t x) {
  assert(ia!=NULL && n>0);
  size_t l=0, r=n-1;
  while(l<r) {
    size_t m = (l+r)/2;
    if(ia[m]<x) l=m+1;
    // else if(ia[m]==x) return m;
    else r=m; // replace with r = m-1 and later return r+1?
  }
  assert(l==r);
  if(ia[l]<x) {
    assert(r==n-1);
    //return n;   // replace with return r+1?
    l = n;
  }
  //\\ printf("binsearch n: %zu x: %ld ris: %zu\n",n,x,l);
  return l;
}



static int uint64_cmp(const void *p, const void *q)
{
  const uint64_t *a = p;
  const uint64_t *b = q;
  
  if(*a < *b) return -1;
  else if(*a > *b) return 1;
  return 0;
}

// read arcs from a text file and store them in a binary file
// removing duplicates
// return the number of arcs saved to the binary file
size_t *arctxt2bin(FILE *in, FILE *out)
{
  assert(in!=NULL && out!=NULL && n!=NULL);
  uint32_t maxarc = 0; // largest arc in the file
  size_t size=10;      // current size of ia[]
  size_t i=0;          // elements in ia[]
  uint64_t *ia = malloc(size*sizeof(*ia));
  if(ia==NULL) quit("arctxt2bin: malloc failed",__LINE__,__FILE__);
    
  int64_t a,b; size_t line=0;  
  while(true) {
    line++;
    int e = fscanf(f,"%" SCNd64 " %" SCNd64,&a,&b);
    if(e==EOF) break;
    // check input
    if(e!=2) {
      fprintf(stderr,"Invalid file content at line %zu",line);
      exit(EXIT_FAILURE);
    }
    if(a<0 || b<0) {
      fprintf(stderr,"Negative arc id at line %zu",line);
      exit(EXIT_FAILURE);
    }
    if(a>UINT32_MAX || b>UINT32_MAX) {
      fprintf(stderr,"Arc id too large at line %zu",line);
      exit(EXIT_FAILURE);
    }
    // update maxarc
    if(a>maxarc) maxarc=a;
    if(b>maxarc) maxarc=b;
    // cobine arcs into a single uint64_t
    uint64_t arc = a<<32 | b;
    // enlarge ia if necessary
    if(i==size) {
        size = size*2;
        ia = realloc(ia,size*sizeof(*ia));
        if(ia==NULL) quit("arctxt2bin: realloc failed",__LINE__,__FILE__);
    }
    assert(size>i);
    ia[i++] = arc;
  }
  // final resize
  size = i;
  ia = realloc(ia,size*sizeof(*ia));
  if(ia==NULL) quit("arctxt2bin: realloc failed",__LINE__,__FILE__);
  fprintf(stderr,"Read %zu arcs\n",size);
  fprintf(stderr,"Max arc id: %u\n",maxarc);
  // sort arcs
  qsort(ia, size, sizeof(*ia), &uint64_cmp);
  // remove duplicates
  size_t j=0; // non duplicate items
  for(i=0;i<size;i++) {
    if(i>0 && ia[i]==ia[i-1]) {
      fprintf(stderr,"Duplicate arc: %ld %ld\n",ia[i]>>32,ia[i]&UINT32_MAX);
      continue;
    }
    ia[j++] = ia[i];
  }
  fprintf(stderr,"Removed %zu duplicates\n",size-j);
  // save to file
  size_t e = fwrite(ia,sizeof(*ia),j,out);
  if(e!=j) quit("arctxt2bin: fwrite failed",__LINE__,__FILE__);
  *n = j;
  return ia;  
}


// -----------------------------

