/* Sum of boolean matrices represented as k2 trees

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
#define default_ext ".sum"

// static functions at the end of the file
static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);



int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname1[PATH_MAX], iname2[PATH_MAX], oname[PATH_MAX];
  #ifdef K2MAT
  char *infofile1=NULL, *infofile2=NULL;
  char *backpfile1=NULL, *backpfile2=NULL; // file with backpointers
  uint32_t rank_block_size = 64; // block size for rank DS  
  #endif
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  bool write = true;
  char *outfile = NULL;
  Use_all_ones_node = true;
  while ((c=getopt(argc, argv, "i:j:r:I:J:o:qhcnvx")) != -1) {
    switch (c) 
      {
      case 'o':
        outfile = optarg; break;                 
      #ifdef K2MAT
      case 'I':
        backpfile1 = optarg; break;                 
      case 'i':
        infofile1 = optarg; break;                 
      case 'J':
        backpfile2 = optarg; break;                 
      case 'j':
        infofile2 = optarg; break;                
      case 'x':
        Use_all_ones_node = false; break;
      case 'r':
        rank_block_size = atoi(optarg); break; // block size of rank structure
      #endif  
      case 'n':
        write = false; break;            
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

  
  // create file names
  sprintf(iname1,"%s",argv[1]);
  sprintf(iname2,"%s",argv[2]);
  if(outfile!=NULL) sprintf(oname,"%s",outfile);
  else       sprintf(oname,"%s%s",argv[1],default_ext); 

  // init matrix variables (valid for b128 and k2t)
  k2mat_t a=K2MAT_INITIALIZER, b=K2MAT_INITIALIZER, ab=K2MAT_INITIALIZER;
  size_t size;

  // load first matrix possibly initializing k2 library
  #ifdef K2MAT
  size = mload_extended(&a, iname1, infofile1, backpfile1, rank_block_size);
  assert(size==a.realsize);
  #else
  size = mload_from_file(&a, iname1);
  #endif
  if (verbose) mshow_stats(&a,iname1,stdout);

  // copy or load second matrix
  size_t size1;
  #ifdef K2MAT
  // possibly load subtree info
  size1 = mload_extended(&b, iname2, infofile2, backpfile2, rank_block_size);
  if(b.fullsize!=a.fullsize) quit("k2 matrices have different internal sizes",__LINE__,__FILE__);
  #else
  size1 = mload_from_file(&b, iname2);
  #endif
  // check sizes correpondds
  if(size1!=size) quit("Input matrices have different sizes",__LINE__,__FILE__);
  if (verbose) mshow_stats(&b,iname2,stdout);

  // do the multiplication show/save the result
  msum(&a,&b,&ab);
  if (verbose || !write) 
    mshow_stats(&ab,oname,stdout);
  if(write) msave_to_file(&ab,oname);

  // free and terminate
  matrix_free(&a);
  matrix_free(&b); // b is distinct from a, deallocate it
  matrix_free(&ab);    
  minimat_reset(); // reset the minimat library and free minimat product table
  // report running time
  fprintf(stderr,"Elapsed time: %.8lf secs\n", ((double) (time(NULL)-start_wc)));
  if(verbose) fprintf(stderr,"==== Done\n");
  return EXIT_SUCCESS;
}


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile1 infile2\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-n         do not write output file, only show stats\n");    
    fprintf(stderr,"\t-o out     outfile name (def. infile1%s)\n",default_ext);
    #ifdef K2MAT
    fprintf(stderr,"\t-i subt1   infile1 subtree info file (subtinfo not used fur sum)\n");
    fprintf(stderr,"\t-j subt2   infile2 subtree info file (subtinfo not used fur sum)\n");
    fprintf(stderr,"\t-I backp1  infile1 backpointers file\n");
    fprintf(stderr,"\t-J backp2  infile2 backpointers file\n");
    fprintf(stderr,"\t-t size    rank block size for k2 compression (def. 64)\n");
    fprintf(stderr,"\t-x         do not compact new 1's submatrices in the result matrix\n");
    #endif  
    fprintf(stderr,"\t-h         show this help message\n");    
    fprintf(stderr,"\t-v         verbose\n\n");
    fprintf(stderr,"Sum two compressed matrices stored in infile1 and infile2\n\n");
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

