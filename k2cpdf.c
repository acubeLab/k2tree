/* Compress a matrix ink2 forat by detecting repeated subtrees 
 * and replacing the root of a repeats subtree with the 0000 value
 * and storing separatly the position of the previous occurrence of the subtree 
 *
 * Since 0000 is used to denot a pruned subtree it cannot be used to
 * mart an submatrix of all 1's  
*/
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
extern bool Use_all_ones_node; // use the special ALL_ONES node in the result
#define default_ext ".ck2"

static void usage_and_exit(char *name);

int main(int argc, char* argv[]) {
  extern char *optarg;
  extern int optind, optopt;
  Use_all_ones_node = false;

  uint32_t rank_block_size = 64;
  uint32_t threshold = 32;
  int verbose = 0, check = 0, write = 1;
  char *infofile1=NULL, *backpfile1=NULL, *outfile=NULL; // file with backpointers

  char* k2name_file = NULL;
  int c;
  while((c = getopt(argc, argv, "hi:I:o:r:t:vcn")) != -1) {
    switch (c) {
      case 'o':
        outfile = optarg; break;                 
      case 'I':
        backpfile1 = optarg; break;                 
      case 'i':
        infofile1 = optarg; break;                 
      case 'r':
        rank_block_size = atoi(optarg); break;
      case 't':
        threshold = atoi(optarg); break;
      case 'h':
        usage_and_exit(argv[0]); break;
      case 'v':
        verbose++; break;
      case 'c':
        check = 1; break;
      case 'n':
        write = 0; break;
      case '?':
        fprintf(stderr, "Unknown option %c\n", optopt);
        exit(1);
    }
  }
  if(verbose > 0) {
    fputs("==== Command line:\n",stdout);
    for(int i=0;i<argc;i++)
      fprintf(stdout," %s",argv[i]);
    fputs("\n",stdout);  
  }

  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;
  k2name_file = argv[1];
  if(k2name_file == NULL) {
    usage_and_exit(argv[0]); 
    exit(1);
  }

  // check command line options
  // threshold >= 4 imples that subtrees of height 1 (root + minimats)
  // are never compressed (recall minimatsize must be == 2) this is a requirement in our algorithms!
  // note that small thredhold can increase the space usage!
  if(threshold % 4 || threshold < 16) {
    fprintf(stderr, "Threshold must be a positive multiple of 4 and at least 16, got %d\n", threshold);
    exit(1);
  }
  if(strlen(k2name_file)> PATH_MAX -16) {
    fprintf(stderr, "Illegal input file name, max size: %d\n", PATH_MAX-16);
    exit(1);
  }

  // create output file name
  char oname[PATH_MAX];
  if(outfile!=NULL) sprintf(oname,"%s",outfile);
  else {
    sprintf(oname,"%s",k2name_file);
    char *dot = strrchr(oname, '.'); // search last dot  
    if (dot && dot != oname) *dot = '\0'; // dot found and it's not the first character, remove ext
    strcat(oname,default_ext);  
  }

  k2mat_t a = K2MAT_INITIALIZER;
  mload_extended(&a, k2name_file, infofile1, backpfile1, rank_block_size); // also init k2 library
  if(minimat_size()!=2) {
    fprintf(stderr, "Error: to compress k2tree mini-matrix size must be 2, got %d\n", minimat_size());
    fprintf(stderr, "(otherwise a spurious 0000 may appear in a mini-matrix)\n");
    matrix_free(&a);
    exit(1);
  }

  // scan input matrix
  size_t pos, nodes, minimats, nz, all1;
  mstats(&a,&pos,&nodes,&minimats,&nz,&all1);

  // check if matrix need to be normalized
  if(a.main_diag_1) 
    puts("## Original matrix has main_diagonal flag, need to be normalized");
  if(a.backp!=NULL) 
    puts("## Original matrix has backpointers, need to be normalized");
  if(all1>0) 
    puts("## Original matrix has ALL_ONES submatrices, need to be normalized");
  if(a.main_diag_1 || a.backp!=NULL || all1>0) {
    puts("## Normalizing input matrix");
    k2mat_t b = mat_zero(&a);
    k2copy_normalise(&a,&b);
    matrix_free(&a);
    a = b;
  }

  // compute number of nonzeros, if verbose show the stats of the original matrix
  if(verbose) {
    printf("Input matrix stats (possibily normalized):\n");
    nz = mshow_stats(&a, basename(k2name_file), stdout);
  }
  else nz = mget_nonzeros(&a);

  // execute subtree compression
  k2mat_t ca = mat_zero(&a);
  k2compress(&a, &ca, threshold, rank_block_size); 

  if(ca.backp==NULL) {
    printf("## No compressible subtree founds. No output file produced\n");
  }
  else {
    if(check || verbose || !write) {
      rank_init(&(ca.r),rank_block_size,&ca); // we need rank info to navigate ca
      size_t nz_ca = mshow_stats(&ca, "Compressed matrix", stdout);
      if(nz_ca!=nz) {
        printf("Compression failed! mismatch in the nuber of nonzeros:\n");
        printf("   %zd (input) vs %zd (compressed)\n",nz,nz_ca);
        exit(1);
      }
      if(check) {
        if(verbose) printf("## Decompressing the k2tree and comparing it with the original\n");
        k2mat_t check_a = mat_zero(&ca);
        size_t pos = 0;
        k2decompress(ca.fullsize, &ca, &pos, &check_a);
        size_t nz_check_a = mshow_stats(&check_a, "Decompressed matrix", stdout);
        if(nz_check_a!=nz) {
          printf("Decompression failed! mismatch in the number of nonzeros:\n");
          printf("   %zd (input) vs %zd (decompressed)\n",nz,nz_check_a);
          exit(1);
        }
        int d = mequals_plain(a.fullsize, &a, &check_a);
        if(d < 0) {
            if(verbose) printf("## Correct decompression!\n");
        } else {
          printf("Error at decompression: generated a different k2tree (error at level: %d)\n",d);
          exit(1);
        }
        matrix_free(&check_a);
      }
    }
    if(write) {
      // save .ck2 file
      msave_to_file(&ca, oname);
      // save .ck2.p file with pointers
      #ifdef SIMPLEBACKPOINTERS
      strcat(oname, ".p");
      #else
      strcat(oname, ".xp");
      #endif
      pointers_write_to_file(ca.backp, oname);
      // NOTE: the rank000 data structure is not stored but recomputed from scratch
    }
  }
  matrix_free(&a);
  matrix_free(&ca);
  minimat_reset();
  return 0;
}

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile.k2\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-o out    outfile name (def. infile1%s)\n",default_ext);
    fprintf(stderr,"\t-I info   infile backpointers file (optional)\n");
    fprintf(stderr,"\t-i info   infile subtree info file (not used)\n");
    fprintf(stderr,"\t-r size   block size for rank 0000 data structure (def. 64)\n");
    fprintf(stderr,"\t-t thrs   smallest subtree size in bits to be removed (def. 32)\n");
    fprintf(stderr,"\t-c        check number of ones in the compressed matrix\n");
    fprintf(stderr,"\t-n        do not write the output file, only show stats\n");
    fprintf(stderr,"\t-h        show this help message\n");    
    fprintf(stderr,"\t-v        verbose\n\n");
    fprintf(stderr,"Compress a k2 tree exploiting the presence of identical subtrees.\n" 
                   "Compute and store in separates files the compressed tree and\n"
                   "its auxiliary information (pointers information)\n\n");
    exit(1);
}
