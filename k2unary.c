/* Demo of unary operations onboolean matrices represented as k2 trees

   Note: internally matrix dimensions are always of the form 2^k times the size 
   of a minimatrix (those stored at the leaves of the tree), with k>0, 
   but the input can be also of a smaller size, and the matrix will be padded with 0's 


   Copyright August 2025-today   ---  giovanni.manzini@unipi.it
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
// definitions to be used for b128 vs k2t-encoded matrices 
#ifdef K2MAT
 #include "k2.h"
 extern bool Use_all_ones_node; // use the special ALL_ONES node in the result
 extern bool Extended_edf;      // compute subtree info on the fly 
#else // definitions for b128 matrices
 #include "b128.h"
 #define K2MAT_INITIALIZER B128MAT_INITIALIZER
 typedef b128mat_t k2mat_t;
 bool Use_all_ones_node; // not used: added for compatibility with k2mat
 bool Extended_edf;      // not used: added for compatibility with k2mat
#endif
// used by both matrix types
#include "bbm.h"
#define default_ext ".unary"

// static functions at the end of the file
static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);


int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname1[PATH_MAX], oname[PATH_MAX];
  #ifdef K2MAT
  char *infofile1=NULL;
  char *backpfile1=NULL; // file with backpointers
  uint32_t rank_block_size = 64; // block size for rank DS  
  #endif
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  char *outfile = NULL;
  Use_all_ones_node = false; Extended_edf = false;
  while ((c=getopt(argc, argv, "i:r:I:o:qhcnv1e")) != -1) {
    switch (c) 
      {
      case 'o':
        outfile = optarg; break;                 
      #ifdef K2MAT
      case 'I':
        backpfile1 = optarg; break;                 
      case 'i':
        infofile1 = optarg; break;                 
      case 'e':
        Extended_edf = true; break; // compute subtree info on the fly                 
      case '1':
        Use_all_ones_node = true; break;
      case 'r':
        rank_block_size = atoi(optarg); break; // block size of rank structure
      #endif          
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
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // always use ALL_ONES root in results
  Use_all_ones_node = true;
  
  // create file names
  sprintf(iname1,"%s",argv[1]);
  if(outfile!=NULL) sprintf(oname,"%s",outfile);
  else       sprintf(oname,"%s%s",argv[1],default_ext); 

  // init matrix variables (valid for b128 and k2tree)
  k2mat_t a=K2MAT_INITIALIZER;
  size_t size, asize;

  // load first matrix possibly initializing k2 library
  #ifdef K2MAT
  size = mload_extended(&asize, &a, iname1, infofile1, backpfile1, rank_block_size);
  #else
  size = mload_from_file(&asize, &a, iname1);
  #endif
  if (verbose) mshow_stats(size,asize,&a,iname1,stdout);

  // add zero matrix 
  k2mat_t b=K2MAT_INITIALIZER,a0=K2MAT_INITIALIZER;
  msum(asize,&b,&a,&a0); // a0 = b+a = 0+a = a
  sprintf(oname,"%s.0.txt",outfile);
  if(verbose)  mshow_stats(size,asize,&a0,oname,stdout);
  mwrite_to_textfile(size,asize, &a0, oname);
  
  // add main diagonal 1's
  madd_identity(&a);
  printf("Caution: madd_identity may add 1's also outside the original matrix size!\n");
  sprintf(oname,"%s.1.txt",outfile);
  mwrite_to_textfile(size,asize, &a, oname);

  // squaring 
  k2mat_t asq = K2MAT_INITIALIZER;
  mmult(asize,&a,&a,&asq);
  if(verbose)   mshow_stats(size, asize,&asq,"(A+I)^2",stdout);
  sprintf(oname,"%s.1sq.txt",outfile);
  msave_to_file(size,asize,&asq,oname);

  // done
  matrix_free(&a);    
  minimat_reset(); // reset the minimat library and free minimat product table
  // report running time
  fprintf(stderr,"Elapsed time: %.8lf secs\n", ((double) (time(NULL)-start_wc)));
  if(verbose) fprintf(stderr,"==== Done\n");
  return EXIT_SUCCESS;
}


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile\n\n",name);
    fprintf(stderr,"Demo of unary operations on the compressed matrices stored in infile\n\n");
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-n        do not write output file, only show stats\n");    
    fprintf(stderr,"\t-o out    outfile name (def. infile1%s)\n",default_ext);
    #ifdef K2MAT
    fprintf(stderr,"\t-1        compact all 1's submatrices in the result matrix\n");
    fprintf(stderr,"\t-i info   infile subtree info file\n");
    fprintf(stderr,"\t-I info   infile backpointers file\n");
    fprintf(stderr,"\t-r size   rank block size for k2 compression (def. 64)\n");
    fprintf(stderr,"\t-e        compute subtree info on the fly (def. no)\n");
    #endif  
    fprintf(stderr,"\t-q        use a single copy when squaring a matrix\n");
    fprintf(stderr,"\t-h        show this help message\n");    
    fprintf(stderr,"\t-v        verbose\n\n");
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

