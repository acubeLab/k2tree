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
  size_t size, asize, totnz=0;
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

  if(threshold % 4) {
    threshold = (threshold + 4 - 1) / 4;
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
  char file_ones[strlen(k2name_file) + 4];
  strcpy(file_ones, k2name_file);
  strcat(file_ones, ".pos");
  mwrite_to_textfile(size, asize, &a, file_ones);

  totnz += mshow_stats(size, asize, &a, basename(k2name_file), stdout);


  k2mat_t ca = K2MAT_INITIALIZER;

  k2compress(asize, &a, &ca, threshold, rank_block); 

  char file_ck2[strlen(k2name_file) + 5];
  strcpy(file_ck2, k2name_file);
  file_ck2[strlen(k2name_file) - 2] = 'c';
  file_ck2[strlen(k2name_file) - 1] = 'k';
  file_ck2[strlen(k2name_file)] = '2';
  file_ck2[strlen(k2name_file) + 1] = '\0';

  if(write) {
    msave_to_file(size, asize, &ca, file_ck2);

    char file_p[strlen(file_ck2) + 4];
    strcpy(file_p, file_ck2);
    strcat(file_p, ".p");
    pointers_write_to_file(ca.p, file_p);

    char file_r[strlen(file_ck2) + 4];
    strcpy(file_r, file_ck2);
    strcat(file_p, ".r");
    rank_write_to_file(ca.r, file_p);
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
        size_t totnz_ca_a = mshow_stats(size, asize, &check_a, "check decompressed tree", stdout);
        if(totnz_ca_a == totnz) {
          int d = mequals(asize, &a, &check_a);
          if(d < 0) {
            if(verbose) printf("Correct decompression!\n");
          } else {
            printf("Error at decompression: generated a different k2tree\n");
          }
        } else {
          printf("Error at decompression: Amount of non zero mismatches! expected %zu, got: %zu\n", totnz, totnz_ca);
        }
      } else {
        printf("Amount of non zero mismatches! expected %zu, got: %zu\n", totnz, totnz_ca);
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
    fprintf(stderr,"\t-b      amount of nodes per block for rank 0000 (def. 64)\n");
    fprintf(stderr,"\t-t      minimum amount of bits to remove a subtree (def. 32)\n");
    fprintf(stderr,"\t-c      check amount of ones of the compressed tree\n");
    fprintf(stderr,"\t-n      do not write the output file, only show stats\n");
    fprintf(stderr,"\t-h      show this help message\n");    
    fprintf(stderr,"\t-v      verbose\n\n");
    fprintf(stderr,"Compute and store in separates files the compressed tree and\n"
                   "its auxiliary information of the input compressed matrix\n"
                   ", based on subtree compression.\n\n");
    exit(1);
}
