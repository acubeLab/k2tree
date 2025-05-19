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

void print_ck2(k2mat_t *a) {
  for(size_t i = 0; i < a->pos; i++) {
    for(size_t bit = 0; bit < 4; bit++) {
      if(i % 2)
        printf("%d", (((a->b[i / 2] >> 4) & (1 << bit)) > 0));
      else
        printf("%d", (((a->b[i / 2] & 15) & (1 << bit)) > 0));
    }
    printf(" ");
  }
  printf("\n");
  if(a->p != NULL) {
    for(size_t i = 0; i < a->p->p_size; i++) {
      printf("%" PRIu32 " ", a->p->p[i]);
    }
    printf("\n");
  }
  if(a->r != NULL) {
    for(size_t i = 0; i < a->r->r_size; i++) {
      printf("%" PRIu32 " ", a->r->r[i]);
    }
  }
  printf("\n");
}

int main(int argc, char* argv[]) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  k2mat_t a = K2MAT_INITIALIZER;
  size_t size, asize, totnz=0;
  uint32_t rank_block = 64;
  uint32_t threshold = 4;

  char* k2name_file = NULL;
  int c;
  while((c = getopt(argc, argv, "b:t:f:")) != -1) {
    switch (c) {
      case 'b':
        rank_block = atoi(optarg); break;
      case 't':
        threshold = atoi(optarg); break;
      case 'f':
        k2name_file = optarg; break;
      case '?':
        fprintf(stderr, "Unknown option %c\n", optopt);
        exit(1);
    }
  }

  if(threshold % 4) {
    fprintf(stderr, "Threshold has to be divisible by 4\n");
    exit(1);
  }

  if(k2name_file == NULL) {
    fprintf(stderr, "Missing file argument: -f <name file>\n");
    exit(1);
  }

  size = mload_from_file(&asize, &a, k2name_file); // also init k2 library
  char file_ones[strlen(k2name_file) + 4];
  strcpy(file_ones, k2name_file);
  strcat(file_ones, ".pos");
  mwrite_to_textfile(size, asize, &a, file_ones);

  totnz += mshow_stats(size, asize, &a, basename(k2name_file), stdout);


  k2mat_t ca = K2MAT_INITIALIZER;

  k2compress(asize, &a, &ca, threshold, rank_block); 

  fprintf(stdout, "COMPRESSED TREE\n");
  fprintf(stdout, "Nodes: %ld\n", ca.pos);

  mshow_stats(size, asize, &ca, basename(k2name_file), stdout);

  char file_ck2[strlen(k2name_file) + 5];
  strcpy(file_ck2, k2name_file);
  file_ck2[strlen(k2name_file) - 2] = 'c';
  file_ck2[strlen(k2name_file) - 1] = 'k';
  file_ck2[strlen(k2name_file)] = '2';
  file_ck2[strlen(k2name_file) + 1] = '\0';
  msave_to_file(size, asize, &ca, file_ck2);

  char file_p[strlen(file_ck2) + 4];
  strcpy(file_p, file_ck2);
  strcat(file_p, ".p");
  pointers_write_to_file(ca.p, file_p);

  char file_r[strlen(file_ck2) + 4];
  strcpy(file_r, file_ck2);
  strcat(file_p, ".r");
  rank_write_to_file(ca.r, file_p);

  matrix_free(&a);
  matrix_free(&ca);
  minimat_reset();
  return 0;
}
