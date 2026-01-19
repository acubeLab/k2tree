/* Compression and decompression of boolean matrices using k2 trees

   Note: internally matrix dimensions are always of the form 2^k times the size 
   of a minimatrix (those stored at the leaves of the tree), with k>0, 
   but the input can be of any size, and the matrix will be padded with 0's 


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
#ifdef K2MAT
#include "k2.h"
#define default_cext ".k2"
extern bool Use_all_ones_node; // use the special ALL_ONES node?
#else
#include "b128.h"
#define default_cext ".b128"
#define K2MAT_INITIALIZER B128MAT_INITIALIZER
typedef b128mat_t k2mat_t;
bool Use_all_ones_node; // not used: added for compatibility with k2mat 
#endif
// used by both matrix type 
#include "bbm.h"
#define default_dext ".bbm"

// static functions at the end of the file
static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);


int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname[PATH_MAX], oname[PATH_MAX];
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  int mmsize = 2;
  bool decompress = false, check = false, write = true;
  char *outfile = NULL;
  Use_all_ones_node = false;
  while ((c=getopt(argc, argv, "o:m:1dchvn")) != -1) {
    switch (c) 
      {
      case 'o':
        outfile = optarg; break;
      case 'm':
        mmsize = atoi(optarg); break;
      case '1':
        Use_all_ones_node = true; break;        
      case 'd':
        decompress = true; break;
      case 'c':
        check = true; break;      
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
  if(check  && decompress)
    quit("Options -c and -d are incompatible",__LINE__,__FILE__);

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // assign default extension and create file names
  char *ext = decompress ? default_dext : default_cext;  
  sprintf(iname,"%s",argv[1]);
  if(outfile!=NULL) sprintf(oname,"%s",outfile);
  else              sprintf(oname,"%s%s",argv[1],ext); 

  k2mat_t a = K2MAT_INITIALIZER;
  size_t size, asize; uint8_t *b = NULL;
  if(decompress) {
    size = mload_from_file(&asize, &a, iname); // also init k2 library
    if (verbose || !write)  
      mshow_stats(size, asize,&a,iname,stdout);
    b= bbm_alloc(size);
    mwrite_to_bbm(b,size,asize, &a);
    if(write) bbm_write(b,size,oname);
  }
  else { // compression
    minimat_init(mmsize);     // init k2 library
    b = bbm_read(iname,&size); // file  ->bbm
    if(verbose>1) bbm_to_ascii(b,size,0,0,size,stdout);
    asize = mread_from_bbm(b,size,&a);
    if (verbose || !write)  
      mshow_stats(size, asize,&a,iname,stdout);
    if(write) msave_to_file(size,asize,&a,oname);  // save k2mat to file
    if(check) {
      uint8_t *bx = bbm_alloc(size);
      mwrite_to_bbm(bx,size,asize, &a);
      ssize_t eq = mequals_bbm(b,size,bx);
      if(eq<0) fprintf(stdout,"Decompressed and input matrix are identical!\n"); 
      else fprintf(stdout,"Decompressed matrix differs at position (%zd,%zd) "
      "dec:%d vs orig:%d\n",eq/size,eq%size, bx[eq], b[eq]);
      free(bx);
    }
  }
  if(verbose>1) bbm_to_ascii(b,size,0,0,size,stdout);
  free(b);
  matrix_free(&a);

  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  if(verbose) fprintf(stderr,"==== Done\n");
  
  return EXIT_SUCCESS;
}

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile\n\n",name);

    fputs("Tool to compress boolean matrices in bbm format (one byte per entry)\n", stderr); 
    #ifdef K2MAT
    fputs("to k2 compressed Plain Depth First format\n\n",stderr);
    #else
    fputs("to B128 (one bit per entry) format\n\n",stderr);
    #endif


    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-d      decompress\n");
    fprintf(stderr,"\t-n      do not write the output file, only show stats\n");
    fprintf(stderr,"\t-o out  outfile name (def. compr: infile%s, decompr: infile%s)\n",
                   default_cext, default_dext);
    #ifdef K2MAT
    fprintf(stderr,"\t-m M    minimatrix size (def. 2) [compression only]\n");
    fprintf(stderr,"\t-1      compact all 1's submatrices [compression only]\n");
    #endif  
    fprintf(stderr,"\t-c      compress->decompress->check\n");
    fprintf(stderr,"\t-h      show this help message\n");    
    fprintf(stderr,"\t-v      verbose\n\n");
    fprintf(stderr,"Default action is to compress filename to filename%s\n\n",default_cext);
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

