/* Compression and decompression of boolean matrices using k2 trees

   k2showinfo: provide info on the encoding of one 
           k2-compressed matrix calling mshow_info 

   For details of the k2 format see how the encoding is done in 
     mread_from textfile in  k2text.c (input is list of nonzero entries)
     mencode_bbm in k2ops.c (input is a one-byte per entry dense matrix)

   This code is ready to be compiled with the B128MAT constant defined
   in that case it will call mshow_info on the b128 compressed matrices

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
#include <unistd.h>
#include <stdbool.h>
#include <libgen.h>
// definitions to be used for b128 vs k2-encoded matrices 
#ifndef K2MAT
#include "b128.h"
#define default_cext ".b128"
#define K2MAT_INITIALIZER B128MAT_INITIALIZER
typedef b128mat_t k2mat_t;
#else // k2mat
#include "k2.h"
#define default_cext ".k2"
#endif

// static functions at the end of the file
static void usage_and_exit(char *name);
// static void quit(const char *msg, int line, char *file);


int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
#ifdef K2MAT
  char *infofile1=NULL, *backpfile1=NULL; // file subtree info and backpointers
  uint32_t rank_block_size = 64; // block size for rank DS  
  #endif
  int c;
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  while ((c=getopt(argc, argv, "hi:I:t:")) != -1) {
    switch (c) 
      {
    #ifdef K2MAT  
      case 'I':
        backpfile1 = optarg; break;                 
      case 'i':
        infofile1 = optarg; break;                 
      case 't':
        rank_block_size = atoi(optarg); break; // block size of rank structure
      #endif
      case 'h':
        usage_and_exit(argv[0]); break;        
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
  }

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  k2mat_t a = K2MAT_INITIALIZER;
  size_t size, asize, totnz=0;
  printf("Matrix file: %s\n",argv[1]);
  #ifdef K2MAT
  size = mload_extended(&a, argv[1], infofile1, backpfile1, rank_block_size);
  #else
  size = mload_from_file(&a, argv[1]); // also init k2 library
  #endif
  asize = a.fullsize;
  fprintf(stdout,"Caution: the following information is incorrect if the input matrix is subtree compressed (ck2 format)\n"); 
  totnz = mshow_stats(&a,basename(argv[1]),stdout);
  puts("");
  matrix_free(&a);
  minimat_reset();

  // statistics
  fprintf(stderr,"Total nonzeros: %zu\n",totnz);
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  
  return EXIT_SUCCESS;
}

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] file1\n\n",name);
    fputs("Options:\n",stderr);
  #ifdef K2MAT
    fprintf(stderr,"\t-i info   subtree info file\n");
    fprintf(stderr,"\t-I info   backpointers file\n");
    fprintf(stderr,"\t-t size   rank block size for k2 compression (def. 64)\n");
    #endif
    fprintf(stderr,"\t-h      show this help message\n");    
    fprintf(stderr,"Get info on the k2-compressed matrix file1\n\n");
    exit(1);
}
