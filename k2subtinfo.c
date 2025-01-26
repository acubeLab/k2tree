/* Computation of subtree information for a k2 matrix 
 *
 * For each node for which the information is provided, we store 
 * the size of the subtrees rooted at its children (except the last child)
 * in order to access any such subtree in constant time. 
 * Since this information is stored into a separate array, we also
 * store the size of such information for any subtree, so that we can 
 * reach in O(1) time such information for any given subtree    
 * 
 * For the details of the enconding, see the long comment before the 
 * function k2dfs_sizes() in k2text.c 
 * 
 * Copyright August 2024-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#undef NDEBUG      // always compile assertions 
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <stdbool.h>
#include <libgen.h>
#include "k2.h"
#define default_ext ".sinfo"
#define BITSxTSIZE 40
#define TSIZEMASK ( (((uint64_t) 1)<<BITSxTSIZE) -1 )


// static functions at the end of the file
static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);
static size_t intsqrt(size_t n);



int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname[PATH_MAX], oname[PATH_MAX];
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  bool check = false, write = true;
  char *outfile = NULL;
  int32_t depth_subtree = 0;
  long node_limit = 0;
  while ((c=getopt(argc, argv, "o:D:N:cnhv")) != -1) {
    switch (c) 
      {
      case 'o':
        outfile = optarg; break;
      case 'c':
        check = true; break;
      case 'n':
        write = false; break;
      case 'D':
        depth_subtree = atoi(optarg); break;
      case 'N':
        node_limit = atol(optarg); break;
      case 'h':
        usage_and_exit(argv[0]); break;        
      case 'v':
        verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
  }
  if(depth_subtree!=0 && node_limit!=0) 
    quit("Options -D and -N are incompatible",__LINE__,__FILE__);
  if(depth_subtree<0 || node_limit<0) 
    quit("Options -D and -N must be non-negative",__LINE__,__FILE__);
    
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

  // assign default extension and create file names
  sprintf(iname,"%s",argv[1]);
  if(outfile!=NULL)sprintf(oname,"%s",outfile);
  else       sprintf(oname,"%s%s",argv[1],default_ext); 

  // read input compressed matrix
  k2mat_t a = K2MAT_INITIALIZER;
  size_t size, asize;
  size = mload_from_file(&asize, &a, iname); // also init k2 library
  if (verbose)  
      mshow_stats(size, asize,&a,iname,stdout);
  // compute encoding information
  vu64_t z;      // resizable array to contain the subtree sizes
  vu64_init(&z);
  uint64_t p;
  size_t pos=0;
  // check if the limit is on the subtree depth or node count
  if(depth_subtree > 0) {
    printf("Computing subtree sizes up to level: %d\n", depth_subtree);
    // visit tree, compute and save subtree sizes in z  
    p = k2dfs_sizes(asize,&a,&pos,&z,depth_subtree);
  }
  else {
    if(node_limit==0) node_limit = intsqrt(a.pos);
    printf("Computing subtree sizes up to node limit: %zu\n", node_limit);
    // visit tree, compute and save subtree sizes in z  
    p = k2dfs_sizes_limit(asize,&a,&pos,&z,(size_t)node_limit);
  }
  assert(pos==a.pos);         // check visit was complete
  assert((p&TSIZEMASK)==a.pos); // low bits contain size of whole matrix 
  printf("Subtree info storage: %zu bytes\n", z.n*sizeof(z.v[0]));
  if(write) {
    FILE *out = fopen(oname,"wb");
    if(!out) quit("Error opening output file",__LINE__,__FILE__);
    vu64_write(out,&z);
    fclose(out);
  }
  if(check || !write) { // if not writing force check 
    size_t znsave = z.n;
    z.n=0; // reset z vector
    pos=0; // start from position 0 in the k2 matrix
    size_t pcheck = k2dfs_check_sizes(asize,&a,&pos,&z,znsave);
    assert(pos==a.pos);      // check visit was complete
    if(z.n!=znsave)
      printf("Subtree storage size mismatch! expected: %zu, got: %zu uint64s\n",znsave,z.n);
    else 
      if(verbose) printf("Subtree storage size matches: %zu uint64s\n",z.n);
    if(a.pos!=pcheck)
      printf("Tree size mismatch! expected: %zu, got: %zu nibbles\n",a.pos,pcheck);
    else 
      if(verbose) printf("Tree size matches: %zu nibbles (%zu bytes)\n",pcheck,(pcheck+1)/2);
  }
  vu64_free(&z);
  matrix_free(&a);
  minimat_reset(); // reset the minimat library and free minimat product table
  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  return EXIT_SUCCESS;
}


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-n      do not write the output file, only show stats and check\n");
    fprintf(stderr,"\t-o out  outfile name (def. infile%s)\n", default_ext);
    fprintf(stderr,"\t-D D    depth limit for storing subtree information (def. do not use depth)\n");
    fprintf(stderr,"\t-N N    #node limit for storing subtree information (def. sqrt(tot_nodes))\n");
    fprintf(stderr,"\t-c      check subtree information\n");
    fprintf(stderr,"\t-h      show this help message\n");    
    fprintf(stderr,"\t-v      verbose\n\n");
    fprintf(stderr,"Compute and store in a separate file the information about the size\n"
                   "of the top subtrees of the input compressed matrix.\n\n");
    exit(1);
}

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

