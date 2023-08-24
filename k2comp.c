/* Compression and decompression of boolean matrices using k2 trees

   Note: internally matrix dimensions are always of the form 2^k times the size 
   of a minimatrix (those stored at the leaves of the tree), with k>0, 
   but the input can be of any size, and the matrix will be padded with 0's 


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
#include "k2.h"
#include "bbm.h"

// default extensions
#define default_cext ".k2"
#define default_dext ".d"

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
  char *ext = NULL;
  while ((c=getopt(argc, argv, "e:m:dcvn")) != -1) {
    switch (c) 
      {
      case 'e':
        ext = optarg; break;        
      case 'm':
        mmsize = atoi(optarg); break;
      case 'd':
        decompress = true; break;
      case 'c':
        check = true; break;      
      case 'n':
        write = false; break;       
      case 'v':
        verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
  }
  if(verbose>0) {
    fputs("==== Command line:\n",stderr);
    for(int i=0;i<argc;i++)
     fprintf(stderr," %s",argv[i]);
    fputs("\n",stderr);  
  }
  if(check  && decompress)
    quit("Options -c and -d are incompatible",__LINE__,__FILE__);
  // assign defualt extension if needed  
  if(ext==NULL) {
    if(decompress) ext = default_dext;
    else ext = default_cext;
  }
  if(strlen(ext)==0)
    quit("Empty extension, cannot overwrite input file",__LINE__,__FILE__);

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // create file names
  sprintf(iname,"%s",argv[1]);
  sprintf(oname,"%s%s",argv[1],ext); 

  int asize;    k2mat_t a = K2MAT_INITIALIZER;
  size_t size; uint8_t *b = NULL;
  if(decompress) {
    size = mload_from_file(&asize, &a, iname); // also init k2 library
    if (verbose || !write)  
      mshow_stats(size, asize,&a,iname,stderr);
    b= bbm_alloc(size);
    mwrite_to_bbm(b,size,asize, &a);
    if(write) bbm_write(b,size,oname);
  }
  else { // compression
    minimat_init(mmsize);     // init k2 library
    b= bbm_read(iname,&size); // file  ->bbm
    if(verbose>1) bbm_to_ascii(b,size,0,0,size,stderr);
    asize = mread_from_bbm(b,size,&a);
    if (verbose || !write)  
      mshow_stats(size, asize,&a,iname,stderr);
    if(write) msave_to_file(size,asize,&a,oname);  // save k2mat to file
    if(check) {
      uint8_t *bx = bbm_alloc(size);
      mwrite_to_bbm(bx,size,asize, &a);
      ssize_t eq = mequals_bbm(b,size,bx);
      if(eq<0) fprintf(stderr,"Decompressed matrix is equal to original!\n"); 
      else fprintf(stderr,"Decompressed matrix differs at position (%zd,%zd) "
      "d:%d vs o:%d\n",eq/size,eq%size, bx[eq], b[eq]);
      free(bx);
    }
  }
  if(verbose>1) bbm_to_ascii(b,size,0,0,size,stderr);
  free(b);
  k2_free(&a);

  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  
  return EXIT_SUCCESS;
}

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] filename\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-c      compress->decompress->check\n");
    fprintf(stderr,"\t-d      decompress\n");
    fprintf(stderr,"\t-n      do not write output file, only show stats\n");
    fprintf(stderr,"\t-m M    minimatrix size (def. 2), compression only\n");
    fprintf(stderr,"\t-e ext  outfile extension (def. compr: \"%s\", decompr: \"%s\")\n",
                   default_cext, default_dext);
    fprintf(stderr,"\t-v      verbose\n\n");
    fprintf(stderr,"Default action is to compress filename to filename.k2x\n\n");
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

