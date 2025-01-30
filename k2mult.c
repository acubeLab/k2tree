/* Multiplication of boolean matrices represented as k2 trees

   Note: internally matrix dimensions are always of the form 2^k times the size 
   of a minimatrix (those stored at the leaves of the tree), with k>0, 
   but the input can be also of a smaller size, and the matrix will be padded with 0's 


   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
// definitions to be used for b128 vs k2-encoded matrices 
#ifdef B128MAT
#include "b128.h"
#define K2MAT_INITIALIZER B128MAT_INITIALIZER
typedef b128mat_t k2mat_t;
bool Use_all_ones_node; // not used: added for compatibility with k2mat
#else // k2mat
#include "k2.h"
extern bool Use_all_ones_node; // use the special ALL_ONES node?
#endif
// used by both matrix type
#include "bbm.h"
#define default_ext ".prod"

// static functions at the end of the file
static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);


int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname1[PATH_MAX], iname2[PATH_MAX], oname[PATH_MAX];
  #ifndef B128MAT
  char *infofile1=NULL, *infofile2=NULL;
  #endif
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  bool check = false, write = true;
  char *outfile = NULL;
  Use_all_ones_node = false;
  bool optimize_squaring = false;    // use a single copy of M to compute M^2
  while ((c=getopt(argc, argv, "i:j:o:qhcnv1")) != -1) {
    switch (c) 
      {
      case 'o':
        outfile = optarg; break;                 
      #ifndef B128MAT
      case 'i':
        infofile1 = optarg; break;                 
      case 'j':
        infofile2 = optarg; break;                 
      case '1':
        Use_all_ones_node = true; break;
      #endif
      case 'c':
        check = true; break;      
      case 'n':
        write = false; break;       
      case 'q':
        optimize_squaring = true; break;       
      case 'h':
        usage_and_exit(argv[0]); break;        
      case 'v':
        verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
  }
  if(verbose>0) {
    fputs("==== Command line:\n",stdout);
    for(int i=0;i<argc;i++)
     fprintf(stdout," %s",argv[i]);
    fputs("\n",stdout);  
  }

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 3) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // check command line options
  if(optimize_squaring && strcmp(argv[1],argv[2])!=0) {
    fprintf(stderr,"Option -q is allowed only when the input matrices are equal\n");
    exit(2);
  }
  #ifndef B128MAT
  if(optimize_squaring && infofile2!=NULL) {
    fprintf(stderr,"Options -q and -j are incompatible (second matrix uses the same info as first)\n");
    exit(3);
  }
  #endif
  
  // create file names
  sprintf(iname1,"%s",argv[1]);
  sprintf(iname2,"%s",argv[2]);
  if(outfile!=NULL) sprintf(oname,"%s",outfile);
  else       sprintf(oname,"%s%s",argv[1],default_ext); 

  // init k2 variables
  k2mat_t a=K2MAT_INITIALIZER, b=K2MAT_INITIALIZER, ab=K2MAT_INITIALIZER;
  size_t size, asize;

  // load first matrix
  size = mload_from_file(&asize, &a, iname1); // also init k2 library
  #ifndef B128MAT
  // possibly load subtree info
  if(infofile1) k2add_subtinfo(&a,infofile1);
  #endif
  if (verbose) mshow_stats(size,asize,&a,iname1,stdout);

  // copy or load second matrix
  if(optimize_squaring) {
    assert(strcmp(iname1,iname2)==0); 
    mmake_pointer(&a,&b); // if b==a, use a pointer and save space
  }
  else {
    size_t bsize, size1 = mload_from_file(&bsize, &b, iname2);
    // possibly load subtree info
    #ifndef B128MAT
    if(infofile2) k2add_subtinfo(&b,infofile2);
    #endif
    // check sizes correpondds
    if(size1!=size) quit("Input matrices have different sizes",__LINE__,__FILE__);
    if(bsize!=asize) quit("k2 matrices have different sizes",__LINE__,__FILE__);
  }
  if (verbose) mshow_stats(size, asize,&b,iname2,stdout);

  // do the multiplication shos/save the result
  mmult(asize,&a,&b,&ab);
  if (verbose || !write) 
    mshow_stats(size, asize,&ab,oname,stdout);
  if(write) msave_to_file(size,asize,&ab,oname);

  // check product if requested: use bbm matrix (n^2 bytes n^3 time) 
  if(check) {
    uint8_t *m2, *m1 = bbm_alloc(size), *m3 = bbm_alloc(size);
    // read m1 
    mwrite_to_bbm(m1,size,asize,&a);
    if(verbose>1) bbm_to_ascii(m1,size,0,0,size,stdout);
    // read m2 if different from m2
    if(strcmp(iname1,iname2)==0) m2=m1;
    else {
      m2 = bbm_alloc(size);
      mwrite_to_bbm(m2,size,asize,&b);
      if(verbose>1) bbm_to_ascii(m2,size,0,0,size,stdout);
    }
    // compute product to m3 = m1 * m2
    // consider using fast_mmult_bbm, but that would require opm
    mmult_bbm(m1,size,m2,m3);
    // uncompress product to m1
    mwrite_to_bbm(m1,size,asize,&ab);
    if(verbose>1) bbm_to_ascii(m1,size,0,0,size,stdout);
    ssize_t eq = mequals_bbm(m1,size,m3);
    if(eq<0) fprintf(stdout,"Product matches the one computed using byte matrices!\n");
    else {
      ssize_t ssize = (ssize_t) size;
      fprintf(stdout,"Product matrix differs at position (%zd,%zd) "
      "prod+uncompr:%d vs uncompr+prod:%d\n",eq/ssize,eq%ssize, m1[eq], m3[eq]);
    }
    free(m1); free(m3);
    if(strcmp(iname1,iname2)) free(m2);
  }

  // free and terminate
  matrix_free(&a);
  if(!optimize_squaring) 
    matrix_free(&b); // b is distinct from a, deallocate it
  matrix_free(&ab);    
  minimat_reset(); // reset the minimat library and free minimat product table
  // report running time
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  return EXIT_SUCCESS;
}


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile1 infile2\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-n        do not write output file, only show stats\n");    
    fprintf(stderr,"\t-o out    outfile name (def. infile1%s)\n",default_ext);
    #ifndef B128MAT
    fprintf(stderr,"\t-1        compact all 1's submatrices in the result matrix\n");
    fprintf(stderr,"\t-i info   infile1 subtree info file\n");
    fprintf(stderr,"\t-j info   infile2 subtree info file\n");
    #endif  
    fprintf(stderr,"\t-q        use a single copy when squaring a matrix\n");
    fprintf(stderr,"\t-c        check multiplication (O(n^3) time and O(n^2) space!)\n");
    fprintf(stderr,"\t-h        show this help message\n");    
    fprintf(stderr,"\t-v        verbose\n\n");
    fprintf(stderr,"Multiply two compressed matrices stored in infile1 and infile2\n\n");
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

