/* Compression and decompression of boolean matrices using k2 trees

   k2sparse: (de)compress matrices in text form: one entry (two indices) per line

   Generates the executable k2sparse.x; if compiled with the B128MAT constant
   defined it generates the executable b128sparse.x which (de)compress matrices 
   in text form to/from the B128 format (one bit x entry)

   For details of the k2 format see how the encoding is done in 
     mread_from_textfile in  k2text.c (input is list of nonzero entries)
     mencode_bbm in k2ops.c (input is a one-byte per entry dense matrix)

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
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
bool Use_all_ones_node; // not used: added for compatibility with k2mat 
#else // defintions for k2 matrices 
#include "k2.h"
#define default_cext ".k2"
extern bool Use_all_ones_node; // use the special ALL_ONES node
#endif
// used by both matrix type 
#define default_dext ".txt"
#define matrix_checker "matrixcmp.x"


// static functions at the end of the file
static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);


int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname[PATH_MAX], oname[PATH_MAX];
  #ifdef K2MAT
  char *infofile1=NULL;
  char *backpfile1=NULL; // file with backpointers
  uint32_t rank_block_size = 64; // block size for rank DS  
  #endif
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  int mmsize = 2;
  bool decompress = false, check = false, write = true;
  int64_t xsize = 0;
  char *outfile = NULL;
  Use_all_ones_node = true;
  while ((c=getopt(argc, argv, "o:m:s:di:I:r:cnhvx")) != -1) {
    switch (c) 
      {
      case 'o':
        outfile = optarg; break;
      case 's':
        xsize = atoll(optarg); break;
      case 'd':
        decompress = true; break;
      case 'c':
        check = true; break;
      case 'n':
        write = false; break;
      #ifdef K2MAT
      case 'm':
        mmsize = atoi(optarg); break;
      case 'I':
        backpfile1 = optarg; break;                 
      case 'i':
        infofile1 = optarg; break;                                 
      case 'x':
        Use_all_ones_node = false; break;
      case 'r':
        rank_block_size = atoi(optarg); break; // block size of rank structure
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

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  if(check  && decompress)
    quit("Options -c and -d are incompatible",__LINE__,__FILE__);
  #ifndef K2MAT
  if(xsize<=0)
    quit("-s parameter is mandatory and must be positive",__LINE__,__FILE__); 
  #else    
  if(xsize<0) 
    quit("-s parameter must be non negative",__LINE__,__FILE__);
  #endif
  // check we are within the hard limit of matrix size
  // note there can be other limitations not tested here
  if(xsize>MaxMatrixSize) 
    #ifndef K2MAT
    quit("Matrix size is too large: see b128.h for the hard limit on size",__LINE__,__FILE__);
    #else    
    quit("Matrix size is too large: see k2.h for the hard limit on size",__LINE__,__FILE__);
    #endif

  // assign default extension and create file names
  char *ext = decompress ? default_dext : default_cext;  
  sprintf(iname,"%s",argv[1]);
  if(outfile!=NULL)sprintf(oname,"%s",outfile);
  else       sprintf(oname,"%s%s",argv[1],ext); 

  k2mat_t a = K2MAT_INITIALIZER;
  size_t size; // size the actual matrix size, asize the internal size 2^k * MMsize
  if(decompress) {
    #ifdef K2MAT
    size = mload_extended(&a, iname, infofile1, backpfile1, rank_block_size);
    #else
    size = mload_from_file(&a, iname);
    #endif
    if (verbose || !write)  
      mshow_stats(&a,iname,stdout);
    if(write) mwrite_to_textfile(&a, oname);
  }
  else { // compression
    minimat_init(mmsize);     // init k2 library
    size_t asize = mread_from_textfile(&a,iname,xsize);
    size = a.realsize;
    assert(asize == a.fullsize);
    assert(xsize==0 || (size==xsize));
    if (verbose || !write)  
      mshow_stats(&a,iname,stdout);
    if(write) msave_to_file(&a,oname);  // save k2mat to file
    if(check) {
      strcat(oname,".check"); // create check file name 
      mwrite_to_textfile(&a, oname); 
      matrix_free(&a);
      // statistics (we are not returning from execlp)
      fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
      if(verbose) fprintf(stderr,"==== Done\n");
      puts("==== Checking compression by calling " matrix_checker);
      char *tmp = strdup(argv[0]);
      char *exedir = dirname(tmp);
      char ename[PATH_MAX];
      sprintf(ename,"%s/%s",exedir, matrix_checker); free(tmp);
      if(verbose>0) fprintf(stderr,"Calling: %s %s %s\n",ename,iname,oname);
      execlp(matrix_checker,ename,iname,oname,NULL);
      quit("Error executiong execlp",__LINE__,__FILE__);
    }
  }
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

    fputs("Tool to (de)compress boolean matrices in sparse text format (one line per entry)\n", stderr); 
    #ifdef K2MAT
    fputs("to k2 compressed Plain Depth First format\n\n",stderr);
    #else
    fputs("to B128 (one bit per entry) format\n\n",stderr);
    #endif

    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-d      decompress to text one line per entry format\n");
    fprintf(stderr,"\t-n      do not write the output file, only show stats\n");
    fprintf(stderr,"\t-o out  outfile name (def. compr: infile%s, decompr: infile%s)\n",
                   default_cext, default_dext);
    #ifndef K2MAT
    fprintf(stderr,"\t-s S    matrix actual size [compression only]\n");
    #else
    fprintf(stderr,"\t-s S    matrix actual size (def. largest index+1) [compression only]\n");
    fprintf(stderr,"\t-m M    minimatrix size (def. 2) [compression only]\n");
    fprintf(stderr,"\t-i info infile subtree info file [decompression only] (def. None)\n");
    fprintf(stderr,"\t-I info infile backpointers file [decompression only] (def. None)\n");
    fprintf(stderr,"\t-r size rank block size for backpointres (def. 64)\n");
    fprintf(stderr,"\t-x      do not compact all 1's submatrices [compression only]\n");
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

