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
#include <math.h>

// definitions to be used for b128 vs k2t-encoded matrices 
#ifdef K2MAT
 #include "k2.h"
 extern bool Use_all_ones_node; // use the special ALL_ONES node in the result
 extern bool Extended_edf;      // compute subtree info on the fly 
 static void quit(const char *msg, int line, char *file);
 static size_t intsqrt(size_t n);
 #else // definitions for b128 matrices
 #include "b128.h"
 #define K2MAT_INITIALIZER B128MAT_INITIALIZER
 typedef b128mat_t k2mat_t;
 bool Use_all_ones_node; // not used: added for compatibility with k2mat
 bool Extended_edf;      // not used: added for compatibility with k2mat
#endif
#define default_ext ".tc"

// static functions at the end of the file
static void usage_and_exit(char *name);

// global from minimats.c
extern uint32_t Minimat_node_ratio;

int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname1[PATH_MAX], oname[PATH_MAX];
  #ifdef K2MAT
  // char *infofile1=NULL;
  // char *backpfile1=NULL; // file with backpointers
  // uint32_t rank_block_size = 64; // block size for rank DS  
  int32_t depth_subtree = 0;
  long node_limit = 0;
  double node_limit_multiplier = 1;
  #endif
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  char *outfile = NULL;
  Use_all_ones_node = true; Extended_edf = false;

  while ((c=getopt(argc, argv, "o:N:M:D:hvxe")) != -1) { 
    // options r:i:I: currently not supported since they would apply to the first iteration only
    switch (c) 
    {
      case 'o':
        outfile = optarg; break;                 
      #ifdef K2MAT
      // case 'I':
      //   backpfile1 = optarg; break;                 
      // case 'i':
      //   infofile1 = optarg; break;                 
      // case 'r':
      //   rank_block_size = atoi(optarg); break; // block size of rank structure
      case 'e':
        Extended_edf = true; break; // compute subtree info on the fly                 
      case 'x':
        Use_all_ones_node = false; break;
      case 'N':
        node_limit = atol(optarg); break;
      case 'M':
        node_limit_multiplier = atof(optarg); break;
      case 'D':  
        depth_subtree = atoi(optarg); break;
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
  #ifdef K2MAT
  if(depth_subtree!=0 && (node_limit!=0 || node_limit_multiplier!=1))  
    quit("Options -D and -N/-M are incompatible",__LINE__,__FILE__);
  if(node_limit!=0 && node_limit_multiplier!=1)  
    quit("Options -N and -M are incompatible",__LINE__,__FILE__);
  if(depth_subtree<0 || node_limit<0) 
    quit("Options -D and -N must be non-negative",__LINE__,__FILE__);
  #endif

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;
  
  // create file names
  sprintf(iname1,"%s",argv[1]);
  if(outfile==NULL) sprintf(oname,"%s.tc.k2",argv[1]);
  else sprintf(oname,"%s",outfile);

  // init matrix variables (valid for b128 and k2tree)
  k2mat_t a=K2MAT_INITIALIZER;
  k2mat_t aIa=a;

  // load first matrix possibly initializing k2 library
  mload_from_file(&a, iname1);
  #ifdef K2MAT
  if(a.main_diag_1) quit("Input matrix with main_diag_1=true not supported (yet)",__LINE__,__FILE__);
  if(a.is_pointer) quit("Input matrix with backpointers not supported (yet)",__LINE__,__FILE__);
  assert(a.subtinfo==NULL);
  assert(a.backp==NULL);
  #endif

  int iter=0;
  while(true) {
    if(verbose) printf("## Iteration %d\n", iter); 
    iter++;
    #ifdef K2MAT
      // compute subtree info and write to file
      vu64_t z;      // resizable array to contain the subtree sizes
      vu64_init(&z);
      uint64_t p;
      size_t pos=0;
      // check if the limit is on the subtree depth or node count
      if(depth_subtree > 0) {
        if(verbose) printf("Computing subtree sizes up to level: %d\n", depth_subtree);
        p = k2dfs_sizes(a.fullsize,&a,&pos,&z,depth_subtree); // visit tree, compute and save subtree sizes in z 
      }
      else {
        if(node_limit==0) node_limit = intsqrt(a.pos);
        node_limit = lround((double)node_limit*node_limit_multiplier); // apply multiplier
        size_t min_node_limit = 1 + 4*Minimat_node_ratio;
        if(node_limit < min_node_limit) node_limit = min_node_limit; // minimum node limit (ensure at least 2 levels)
        if(verbose) printf("Computing sizes for subtrees larger than %zu nodes\n", node_limit);
        p = k2dfs_sizes_limit(a.fullsize,&a,&pos,&z,(size_t)node_limit); // visit tree, compute and save subtree sizes in z 
        node_limit =0; // reset for next iterqation
      }
      a.subtinfo_size = z.n; 
      a.subtinfoarray = a.subtinfo = z.v; // save subtree info in a
      assert(pos==a.pos);         // check visit was complete
      assert((p&TSIZEMASK)==a.pos); // low bits contain size of whole matrix 
    #endif  
    if(verbose) mshow_stats(&a,"Current matrix",stdout);
    #ifdef K2MAT
    k2mat_t b=K2MAT_INITIALIZER;
    if(verbose) printf("Adding identity\n");
    mmake_pointer(&a,&b); 
    madd_identity(&b); // b = a + I
    if(verbose) printf("Multiplying (A+I)A\n");
    mmult(&b,&a,&aIa); // aIa = (a+I)a
    #else
    if(verbose) printf("Multiplying A*A\n");
    mmult(&a,&a,&aIa); // aIa  now is A^2
    if(verbose) printf("Computing A*A + A\n");
    madd(&aIa,&a); // aIa = A^2 + A
    #endif    
    if(verbose) printf("Checking if fixed point\n");
    if(mequals(&a,&aIa)==true) break;
    k2mat_t temp = a; a = aIa; aIa = temp; // swap a and aIa for the next iteration
    matrix_free(&aIa);   
  }
  if(verbose) printf("Fixed point reached, stopping\n");
  if(verbose)  mshow_stats(&aIa,"Transitive closure matrix",stdout);  
  msave_to_file(&aIa, oname);
  // done
  matrix_free(&aIa);  
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
    fputs("Options:\n",stderr); 
    fprintf(stderr,"\t-o out    outfile name (def. infile%s)\n",default_ext);
    #ifdef K2MAT
    fprintf(stderr,"\t-e        compute subtree info on the fly (def. no)\n");
    fprintf(stderr,"\t-x        do not compact new 1's submatrices in the result matrix\n");    
    fprintf(stderr,"\t-D D      depth limit for subtree information (def. ignore depth)\n");
    fprintf(stderr,"\t-N N      node limit for subtree information (def. sqrt(tot_nodes))\n");
    fprintf(stderr,"\t-M M      multiplier for node limit (def. 1)\n");
    #endif  
    fprintf(stderr,"\t-h        show this help message\n");    
    fprintf(stderr,"\t-v        verbose\n\n");
    fprintf(stderr,"Compute the transitive of the compressed matrix stored in infile\n\n");
    exit(1);
}

#ifdef K2MAT
// compute integer square root
static size_t intsqrt(size_t n) {
  // assert(n>=0);
  size_t x = n;
  size_t y = (x + 1) / 2;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2;
  }
  assert(x*x <= n && (x+1)*(x+1) > n);
  return x;
}

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
#endif

