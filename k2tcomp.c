/* Compression and decompression of boolean matrices using k2 trees

   k2tcomp: (de)compress matrices in text form: one entry (two indices) per line

   This code is ready to be compiled with the B128MAT constant defined
   in that case it will (de)encode the matrices in text form to/from 
   the B128 format, which is a dense format with one bit per entry. 
   This option is currently noty supported

   For the details of the k2 format compression see k2text.c

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
#ifdef B128MAT
#include "b128.h"
#define default_cext ".b128"
#define K2MAT_INITIALIZER B128MAT_INITIALIZER
typedef b128mat_t k2mat_t;
#else // k2mat
#include "k2.h"
#define default_cext ".k2"
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
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  int mmsize = 2;
  bool decompress = false, check = false, write = true;
  int64_t xsize = 0;
  char *ext = NULL;
  while ((c=getopt(argc, argv, "e:m:s:dcnhv")) != -1) {
    switch (c) 
      {
      case 'e':
        ext = optarg; break;
      case 'm':
        mmsize = atoi(optarg); break;
      case 's':
        xsize = atoll(optarg); break;
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
  // assign default extension if needed  
  if(ext==NULL) {
    if(decompress) ext = default_dext;
    else ext = default_cext;
  }
  if(strlen(ext)==0)
    quit("Empty extension, cannot overwrite input file",__LINE__,__FILE__);
  #ifdef B128MAT
  if(xsize<=0)
    quit("-s parameter is mandatory and must be positive",__LINE__,__FILE__); 
  #else    
  if(xsize<0) 
    quit("-s parameter must be non negative",__LINE__,__FILE__);
  #endif

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // create file names
  sprintf(iname,"%s",argv[1]);
  sprintf(oname,"%s%s",argv[1],ext); 

  int asize;    k2mat_t a = K2MAT_INITIALIZER;
  size_t size, asizetmp;
  if(decompress) {
    size = mload_from_file(&asizetmp, &a, iname); // also init k2 library
    asize = asizetmp; // ugly hack untill we switch every size to size_t
    if (verbose || !write)  
      mshow_stats(size, asize,&a,iname,stdout);
    if(write) mwrite_to_textfile(size,asize, &a, oname);
  }
  else { // compression
    minimat_init(mmsize);     // init k2 library
    size_t asize = mread_from_textfile(&size,&a,iname,xsize);
    assert(xsize==0 || (size==xsize));
    if (verbose || !write)  
      mshow_stats(size, asize,&a,iname,stdout);
    if(write) msave_to_file(size,asize,&a,oname);  // save k2mat to file
    if(check) {
      sprintf(oname,"%s%s.check",argv[1],ext); // create check file name  
      mwrite_to_textfile(size,asize, &a, oname);
      matrix_free(&a);
      puts("==== Checking compression by calling " matrix_checker);
      char *tmp = strdup(argv[0]);
      char *exedir = dirname(tmp);
      char ename[PATH_MAX];
      sprintf(ename,"%s/%s",exedir, matrix_checker); free(tmp);
      if(verbose>0) fprintf(stderr,"Calling: %s %s %s\n",ename,iname,iname);
      execlp(matrix_checker,ename,iname,oname,NULL);
      quit("Error executiong execlp",__LINE__,__FILE__);
    }
  }
  matrix_free(&a);

  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  
  return EXIT_SUCCESS;
}

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] filename\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-d      decompress\n");
    fprintf(stderr,"\t-n      do not write the output file, only show stats\n");
    #ifdef B128MAT
    fprintf(stderr,"\t-s S    matrix actual size, compression only\n");
    #else
    fprintf(stderr,"\t-s S    matrix actual size (def. largest index), compression only\n");
    #endif  
    fprintf(stderr,"\t-m M    minimatrix size (def. 2), compression only\n");
    fprintf(stderr,"\t-e ext  outfile extension (def. compr: \"%s\", decompr: \"%s\")\n",
                   default_cext, default_dext);
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

