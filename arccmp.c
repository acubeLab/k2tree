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
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <time.h>


// prototypes of static functions
static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);

// global verbosity level
int Verbose=0;


// in a sorted uint64_t array ia[0,n-1] containing distinct values find 
// the first entry >= x using binary search 
// return n if no such entry exists
static size_t binsearch(uint64_t *ia, size_t n, uint64_t x) {
  assert(ia!=NULL && n>0);
  size_t l=0, r=n-1;
  while(l<r) {
    size_t m = (l+r)/2;
    if(ia[m]<x) l=m+1;
    //\\ else if(ia[m]==x) return m;
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
uint64_t *arctxt2bin(FILE *f, size_t *n)
{
  assert(f!=NULL && n!=NULL);
  uint32_t maxarc = 0; // largest arc in the file
  size_t size=10;      // current size of ia[]
  size_t i=0;          // elements in ia[]
  uint64_t *ia = malloc(size*sizeof(*ia));
  if(ia==NULL) quit("arctxt2bin: malloc failed",__LINE__,__FILE__);
    
  int64_t a,b; size_t line=0;  
  while(true) {
    int e = fscanf(f,"%" SCNd64 " %" SCNd64,&a,&b);
    if(e==EOF) break;
    line++;
    // check input
    if(e!=2) {
      fprintf(stderr,"Invalid file content at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(a<0 || b<0) {
      fprintf(stderr,"Negative arc id at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(a>UINT32_MAX || b>UINT32_MAX) {
      fprintf(stderr,"Arc id too large at line %zu\n",line);
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
  if(Verbose>0) {fprintf(stderr,"< Read %zu arcs\n",size);
                 fprintf(stderr,"< Max arc id: %u\n",maxarc);}
  // sort arcs
  qsort(ia, size, sizeof(*ia), &uint64_cmp);
  // remove duplicates
  size_t j=0; // non duplicate items
  for(i=0;i<size;i++) {
    if(i>0 && ia[i]==ia[i-1]) {
      if(Verbose>0) fprintf(stderr,"< Duplicate arc: %ld %ld\n",ia[i]>>32,ia[i]&UINT32_MAX);
      continue;
    }
    ia[j++] = ia[i];
  }
  if(Verbose>0) {fprintf(stderr,"< Removed %zu duplicates\n",size-j);
                 fprintf(stderr,"< %zu arcs remaining\n",j);}
  // return
  *n = j;
  return ia;  
}


size_t arctxtcheck(FILE *f,uint64_t m1[], size_t n) {
  assert(f!=NULL && m1!=NULL);
  assert(n>0);
  uint32_t maxarc = 0; // largest arc in the file
  size_t err = 0;
  // allocate and clean bit array 
  uint64_t *bits = calloc((n+63)/64,8);
  if(bits==NULL) quit("arctxtcheck: calloc failed",__LINE__,__FILE__); 

  // main loop reading the matrix 2 file
  int64_t a,b; size_t line=0, dup=0;  
  while(true) {
    int e = fscanf(f,"%" SCNd64 " %" SCNd64,&a,&b);
    if(e==EOF) break;
    line++;
    // check input
    if(e!=2) {
      fprintf(stderr,"Invalid file content at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(a<0 || b<0) {
      fprintf(stderr,"Negative arc id at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(a>UINT32_MAX || b>UINT32_MAX) {
      fprintf(stderr,"Arc id too large at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    // update maxarc
    if(a>maxarc) maxarc=a;
    if(b>maxarc) maxarc=b;
    // combine arcs into a single uint64_t
    uint64_t arc = a<<32 | b;

    size_t pos = binsearch(m1,n,arc);
    if(pos==n || m1[pos]!=arc) {
      fprintf(stderr,"> unmatched %ld %ld\n",a,b);
      err++;
    }
    else { // check for matching duplicates in m2
      assert(m1[pos]==arc);
      if(bits[pos/64]&(1UL<<(pos%64))) {
        if(Verbose>0) fprintf(stderr,"> Duplicate arc: %ld %ld\n",a,b);
        dup++;
        continue;
      }
      bits[pos/64] |= (1UL<<(pos%64));
      assert( (bits[pos/64]&(1UL<<(pos%64))) );
    }
  }
  if(Verbose>0) {
    fprintf(stderr,"> Read %zu arcs\n",line);
    fprintf(stderr,"> Found %zu duplicates\n",dup);
    fprintf(stderr,"> Max arc id: %u\n",maxarc);
  }
  for(size_t i=0;i<n;i++)
    if( (bits[i/64]&(1UL<<(i%64))) ==0) {
      fprintf(stderr,"< unmatched %ld %ld\n",m1[i]>>32,m1[i]&UINT32_MAX);
      err++;
    }
  free(bits);
  return err;  
}

// -----------------------------



int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int c;
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  while ((c=getopt(argc, argv, "hv")) != -1) {
    switch (c) 
      {
      case 'h':
        usage_and_exit(argv[0]); break;        
      case 'v':
        Verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
  }
  if(Verbose>0) {
    fputs("==== Command line:\n",stdout);
    for(int i=0;i<argc;i++)
     fprintf(stderr," %s",argv[i]);
    fputs("\n",stderr);  
  }

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 3) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // open files
  FILE *f1 = fopen(argv[1],"rt");
  if(f1==NULL) quit("Unable to open matrix file 1",__LINE__,__FILE__);
  FILE *f2 = fopen(argv[2],"rt");
  if(f2==NULL) quit("Unable to open matrix file 2",__LINE__,__FILE__);

  fprintf(stderr,"Reading matrix 1\n");
  size_t n;
  uint64_t *m1 = arctxt2bin(f1,&n);
  fclose(f1);
  if(n==0) {
    fprintf(stderr,"Matrix 1 has no arcs\n");
    free(m1);
    fclose(f2);
    return EXIT_FAILURE;
  }
  fprintf(stderr,"Reading matrix 2\n");
  size_t err = arctxtcheck(f2,m1,n);
  free(m1);
  fclose(f2);
  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  // exit
  if(err>0) {
    printf("%zu mismatches found\n",err);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] matrix1 matrix2\n\n",name);
    fputs("Compare arcs in text files matrix1 and matrix2\n",stderr);
    fputs("Reporting mismatches and duplicates\n",stderr);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-h      show this help message\n");
    fprintf(stderr,"\t-v      verbose\n\n");
    exit(1);
}


// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

