/* Computation of subtree information for a k2 matrix 
 *
 * For each node for which the information is provided, we store 
 * the size of the subtrees rooted at its children (except the last child)
 * in order to access any such subtree in constant time. 
 * Since this information is stored into a separate array, we also
 * store the size of such information for any subtree, so that we can 
 * reach in O(1) time such information for any given subtree    
 * 
 * For the details of the encoding, see the long comment before the 
 * function k2dfs_sizes() in k2text.c [change this]
 * 
 * Extended version (SIMPLEBACKPOINTERS not defined, currently not used):
 * If the option -p is used, we assume the nibble 0000 denotes a 
 * pointer to a subtree. The destination of the pointer subtree is stored
 * in the array of uint64_t m->backp which is read from pinfile. 
 * The destination of the i-th (in left-to-right order) pointer is 
 * stored in the lower BITSxTSIZE (40) bits
 * of backp->node[i]. After computing the subtree information we do a second
 * pass where we store in the remaining (24) bits the starting position of 
 * the subtree information for that subtree (if any information is available).
 * This use of backpointers requires that for each subtree it is available 
 * the size of the subtree for all children (not excluding the last one)
 * The subtree info file created by this version has default extension .xsinfo
 * and the backpointer info file has default extension .xpinfo (use -P to change it)
 * 
 * Version with SIMPLEBACKPOINTERS defined:
 * We do not store the subtree information for the destination of backpointers
 * (the rationale is that this would be availble only for large repeated submatrices). 
 * As a consequence:
 *   1. we use an uint32_t for each backpointer (see pointers.h)
 *   2. we do not store the subtree information for the last child
 * The subtree info file has default extension .sinfo
 * This version uses less space, but it is less efficient since
 * when accessing a subtree which is destination of a backpointer
 * we do not have the subtree information.
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
#include <math.h>
#include "k2.h"
#ifdef SIMPLEBACKPOINTERS
#define default_ext ".sinfo"     // extension for simple (old style) subtree info
#define default_pext ".error"    // no such file for simple backpointers
#else
#define default_ext ".xsinfo"
#define default_pext ".xpinfo"
#endif

#ifndef SIMPLEBACKPOINTERS
#pragma message "Compiling with FULL 64bit BACKPOINTERS and subtree information for all children"
#endif

// global from minimats.c
extern uint32_t Minimat_node_ratio;

// static functions at the end of the file
static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);
static size_t intsqrt(size_t n);

int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname[PATH_MAX], oname[PATH_MAX], poname[PATH_MAX];
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  bool check = false, write = true;
  char *outfile=NULL, *pinfile = NULL, *poutfile=NULL;
  int32_t depth_subtree = 0;
  long node_limit = 0;
  double node_limit_multiplier = 1;

  #ifdef SIMPLEBACKPOINTERS
  while ((c=getopt(argc, argv, "o:D:N:M:cnhv")) != -1) {
  #else
  while ((c=getopt(argc, argv, "o:p:P:N:M:cnhv")) != -1) {
  #endif
    switch (c) 
      {
      case 'o':
        outfile = optarg; break;
      case 'P':
        poutfile = optarg; break;        
      case 'p':
        pinfile = optarg; break;
      case 'c':
        check = true; break;
      case 'n':
        write = false; break;
      case 'D':  
        depth_subtree = atoi(optarg); break;
      case 'N':
        node_limit = atol(optarg); break;
      case 'M':
        node_limit_multiplier = atof(optarg); break;
      case 'h':
        usage_and_exit(argv[0]); break;        
      case 'v':
        verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
  }
  if(depth_subtree!=0 && (node_limit!=0 || node_limit_multiplier!=1))  
    quit("Options -D and -N/-M are incompatible",__LINE__,__FILE__);
  if(node_limit!=0 && node_limit_multiplier!=1)  
    quit("Options -N and -M are incompatible",__LINE__,__FILE__);
  if(depth_subtree<0 || node_limit<0) 
    quit("Options -D and -N must be non-negative",__LINE__,__FILE__);
  if(poutfile!=NULL && pinfile==NULL) 
    quit("Option -P requires option -p",__LINE__,__FILE__);

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
  if(outfile!=NULL) sprintf(oname,"%s",outfile);
  else       sprintf(oname,"%s%s",argv[1],default_ext); 
  if(poutfile!=NULL) sprintf(poname,"%s",poutfile);
  else sprintf(poname,"%s%s",argv[1],default_pext);
 
  // read input matrix (could be subtree compressed, we cannot know...)
  k2mat_t a = K2MAT_INITIALIZER;
  mload_from_file(&a, iname); // also init k2 library
  // show information acquired so far from the input files 
  if (verbose) {
    if(a.is_pointer)
      fprintf(stdout,"Caution: the stats information is incorrect: the matrix is subtree compressed\n"); 
    mshow_stats(&a,iname,stdout);
  }
  else if(a.is_pointer)
     fprintf(stdout,"Caution: the matrix is subtree compressed (not relevant here)\n");
  a.is_pointer =false; // reset flag

  // compute subtree information
  vu64_t z;      // resizable array to contain the subtree sizes
  vu64_init(&z);
  uint64_t p;
  size_t pos=0;
  // check if the limit is on the subtree depth or node count
  if(depth_subtree > 0) {
    if(verbose) printf("Computing subtree sizes up to level: %d\n", depth_subtree);
    // visit tree, compute and save subtree sizes in z  
    p = k2dfs_sizes(a.fullsize,&a,&pos,&z,depth_subtree);
  }
  else {
    if(node_limit==0) node_limit = intsqrt(a.pos);
    node_limit = lround((double)node_limit*node_limit_multiplier); // apply multiplier
    size_t min_node_limit = 1 + 4*Minimat_node_ratio;
    if(node_limit < min_node_limit) node_limit = min_node_limit; // minimum node limit (ensure at least 2 levels)
    if(verbose) printf("Computing sizes for subtrees larger than %zu nodes\n", node_limit);
    // visit tree, compute and save subtree sizes in z  
    p = k2dfs_sizes_limit(a.fullsize,&a,&pos,&z,(size_t)node_limit);
  }
  assert(pos==a.pos);         // check visit was complete
  assert((p&TSIZEMASK)==a.pos); // low bits contain size of whole matrix 
  if(verbose) printf("Subtree info takes: %zu bytes\n", z.n*sizeof(z.v[0]));

  if(check || !write) { // if not writing force check of subtree information
    size_t znsave = z.n;
    z.n=0; // reset z vector
    pos=0; // start from position 0 in the k2 matrix
    size_t pcheck = k2dfs_check_sizes(a.fullsize,&a,&pos,&z,znsave);
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


  #ifndef SIMPLEBACKPOINTERS
  // enrich backpointers if option -P was used to provide them 
  // this implies that the matrix was compressed with 64-bit pointers
  // and we add in the high 24 bits the starting position of the subtree info
  // see the comment "Extended version" at the beginning of this file
  // Note rank information is not needed since we do not follow pointers
  if(pinfile!=NULL) {
    a.backp = pointers_load_from_file(pinfile);
    if(a.backp==NULL) quit("Error loading pointer information",__LINE__,__FILE__);
    // compute stored order of the pointers, and init index to 0
    pointers_sort(a.backp);
    // compute info
    size_t znsave = z.n;
    z.n=0; // reset z vector
    pos=0; // start from position 0 in the k2 matrix
    k2dfs_compute_backpointer_info(a.fullsize,&a,&pos,&z);
    assert(pos==a.pos);      // check visit was complete
    assert(z.n==znsave);     // check that we have scanned z
  }
  #endif

  // if -n was not used save the computed structures 
  if(write) {
    FILE *out = fopen(oname,"wb");
    if(!out) quit("Error opening output file",__LINE__,__FILE__);
    vu64_write(out,&z);
    fclose(out);
    if(verbose) fprintf(stdout,"Subtree information written to file: %s\n",oname);
    if(pinfile!=NULL) {
      pointers_write_to_file(a.backp,poname); // write pointer information if available
      if(verbose) fprintf(stdout,"Enriched pointer information written to file: %s\n",poname);
    }
  }
  // free resources
  vu64_free(&z);
  matrix_free(&a);
  minimat_reset(); // reset the minimat library and free minimat product table
  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  if(verbose) fprintf(stderr,"==== Done\n");
  return EXIT_SUCCESS;
}


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-o out  outfile name (def. infile%s)\n", default_ext);
    fprintf(stderr,"\t-n      do not write the output file(s), only show stats and check\n");
    fprintf(stderr,"\t-D D    depth limit for subtree information (def. ignore depth)\n");
    fprintf(stderr,"\t-N N    node limit for subtree information (def. sqrt(tot_nodes))\n");
    fprintf(stderr,"\t-M M    multiplier for node limit (def. 1)\n");
    #ifdef SIMPLEBACKPOINTERS
    fprintf(stderr,"\t-D D    depth limit for subtree information (def. ignore depth)\n");
    #else
    fprintf(stderr,"\t-p pin  file containing backpointer information\n");
    fprintf(stderr,"\t-P pout outfile for backpointer-subtree information (def. infile%s)\n", default_pext);
    #endif
    fprintf(stderr,"\t-c      check subtree information\n");
    fprintf(stderr,"\t-h      show this help message\n");    
    fprintf(stderr,"\t-v      verbose\n\n");
    fprintf(stderr,"Compute and store in a separate file the information about the size\n"
                   "of the largest subtrees of the input matrix.\n"
                   #ifdef SIMPLEBACKPOINTERS
                    "\n");
                   #else
                   "If -p option is used backpointers are enriched with subtree information.\n\n");
                   #endif
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

