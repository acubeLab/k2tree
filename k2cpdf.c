// c includes
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <string.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <stdbool.h>
#include <libgen.h>
#include <inttypes.h>
#include <stdio.h> 

// local includes
#include "k2.h"

static void usage_and_exit(char *name);

int main(int argc, char* argv[]) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  k2mat_t a = K2MAT_INITIALIZER;
  size_t size, asize;
  uint32_t rank_block = 64;
  uint32_t threshold = 32;
  int verbose = 0, check = 0, write = 1;

  char* k2name_file = NULL;
  int c;
  while((c = getopt(argc, argv, "b:t:hvcn")) != -1) {
    switch (c) {
      case 'b':
        rank_block = atoi(optarg); break;
      case 't':
        threshold = atoi(optarg); break;
      case 'h':
        usage_and_exit(argv[0]); break;
      case 'v':
        verbose = 1; break;
      case 'c':
        check = 1; break;
      case 'n':
        write = 0; break;
      case '?':
        fprintf(stderr, "Unknown option %c\n", optopt);
        exit(1);
    }
  }

  // threshold >= 4 imples that subtrees of height 1 (root + minimats)
  // are never compressed (recall minimatsize must be == 2)
  if(threshold % 4 || threshold < 4) {
    fprintf(stderr, "Threshold must be a positive multiple of 4 and at least 16, got %d\n", threshold);
    exit(1);
  }

  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;
  k2name_file = argv[1];
  if(k2name_file == NULL) {
    usage_and_exit(argv[0]); 
    exit(1);
  }

  if(verbose > 0) {
    fputs("==== Command line:\n",stdout);
    for(int i=0;i<argc;i++)
      fprintf(stdout," %s",argv[i]);
    fputs("\n",stdout);  
  }

  size = mload_from_file(&asize, &a, k2name_file); // also init k2 library
  if(minimat_size()!=2) {
    fprintf(stderr, "Error: to compress k2tree mini-matrix size must be 2, got %d\n", minimat_size());
    fprintf(stderr, "(otherwise a spurious 0000 may appear in a mini-matrix)\n");
    matrix_free(&a);
    exit(1);
  }

  // compute number of nonzeros, if verbose show the stats of the original matrix
  size_t totnz = 0;
  if(verbose) {
    printf("Original matrix stats:\n");
    totnz = mshow_stats(size, asize, &a, basename(k2name_file), stdout);
  }
  else totnz = mget_nonzeros(asize, &a);

  // execute subtree compression
  k2mat_t ca = K2MAT_INITIALIZER;
  k2compress(asize, &a, &ca, threshold, rank_block); 


  if(ca.backp==NULL) {
    printf("No compressible subtree founds. No output file produced\n");
  }
  else {

    //!!! TODO: fix this, not nice 
    // replace last two chars of input file name (likely k2) with ck2
    char file_ck2[strlen(k2name_file) + 5];
    strcpy(file_ck2, k2name_file);
    file_ck2[strlen(k2name_file) - 2] = 'c';
    file_ck2[strlen(k2name_file) - 1] = 'k';
    file_ck2[strlen(k2name_file)] = '2';
    file_ck2[strlen(k2name_file) + 1] = '\0';

    if(write) {
      // save .ck2 file
      msave_to_file(size, asize, &ca, file_ck2);
      // save .ck2.p file with pointers
      char file_p[strlen(file_ck2) + 5];
      strcpy(file_p, file_ck2);
      #ifdef SIMPLEBACKPOINTERS
      strcat(file_p, ".p");
      #else
      strcat(file_p, ".xp");
      #endif
      pointers_write_to_file(ca.backp, file_p);

      // NOTE: the rank000 data structure is not stored but recomputed from scratch
    }
    if(check || verbose || !write) {
      size_t totnz_ca = 0;
      totnz_ca = mshow_stats(size, asize, &ca, basename(file_ck2), stdout);
      if(check) {
        if(totnz_ca == totnz) {
          k2mat_t check_a = K2MAT_INITIALIZER;
          size_t pos = 0;
          if(verbose) printf("Decompressing the k2tree and comparing it with the original\n");
          k2decompress(asize, &ca, &pos, &check_a);
          size_t totnz_ca_a = mshow_stats(size, asize, &check_a, "Decompressed matrix", stdout);
          if(totnz_ca_a == totnz) {
            int d = mequals(asize, &a, &check_a);
            if(d < 0) {
              if(verbose) printf("Correct decompression!\n");
            } else {
              printf("Error at decompression: generated a different k2tree (error at level: %d)\n",d);
            }
          } else {
            printf("Error at decompression: Amount of non zero mismatches! expected %zu, got: %zu\n", totnz, totnz_ca);
          }
        } else {
          printf("Number of nonzeros mismatches! expected %zu, got: %zu\n", totnz, totnz_ca);
        }
      }
    }
  }

  matrix_free(&a);
  matrix_free(&ca);
  minimat_reset();
  return 0;
}

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-b      block size for rank 0000 data structure (def. 64)\n");
    fprintf(stderr,"\t-t      smallest subtree size in bits to be removed (def. 32)\n");
    fprintf(stderr,"\t-c      check number of ones in the compressed matrix\n");
    fprintf(stderr,"\t-n      do not write the output file, only show stats\n");
    fprintf(stderr,"\t-h      show this help message\n");    
    fprintf(stderr,"\t-v      verbose\n\n");
    fprintf(stderr,"Compress a k2 tree exploiting the presence of identical subtrees.\n" 
                   "Compute and store in separates files the compressed tree and\n"
                   "its auxiliary information (pointers information)\n\n");
    exit(1);
}
