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
#define InputFile_ext ".bin"
#define OutFile_ext ".k2x"


static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);





void show_stats(int size, k2mat_t *a, char *name) {
  size_t pos, nodes, minimats;
  int levels = mstats(size,a,&pos,&nodes,&minimats);
  printf("%s -- Levels: %d, Pos: %zd, Nodes: %zd, Minimats: %zd\n",
         name,levels,pos,nodes,minimats);
}



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
  bool decompress = false, check = false;
  char *iext = InputFile_ext;
  char *oext = OutFile_ext;
  while ((c=getopt(argc, argv, "i:o:m:dcv")) != -1) {
    switch (c) 
      {
      case 'i':
        iext = optarg; break;        
      case 'o':
        oext = optarg; break;        
      case 'm':
        mmsize = atoi(optarg); break;
      case 'd':
        decompress = true; break;
      case 'c':
        check = true; break;      
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


  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // create file names
  sprintf(iname,"%s%s",argv[1],iext);
  sprintf(oname,"%s%s",argv[1],oext); 

  // init k2 library
  minimat_init(mmsize);
  k2mat_t a = K2MAT_INITIALIZER;
  int size, asize = 0;
  uint8_t *b = NULL;
  if(decompress) {
    size = mload(&asize, &a, iname);
    if (verbose) show_stats(asize,&a,"A");
    b= bbm_alloc(size);
    mwrite_to_bbm(b,size,asize, &a);
    bbm_write(b,size,oname);
  }
  else {
    b= bbm_read(iname,&size); // file  ->bbm
    if(verbose>1) bbm_to_ascii(b,size,0,0,size,stderr);
    asize = mread_from_bbm(b,size,&a);
    if (verbose) show_stats(asize,&a,"A");
    msave(size,asize,&a,oname);
    if(check) {
      uint8_t *bx = bbm_alloc(size);
      mwrite_to_bbm(bx,size,asize, &a);
      bool eq = mequals_bbm(b,size,bx);
      free(bx);
      if(eq) fprintf(stderr,"Decompressed matrix is equal to original!\n"); 
      else  quit("Decompressed matrix is different from original",__LINE__,__FILE__);
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
    fprintf(stderr,"Usage:\n\t  %s [options] basename\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-c      compress->decompress->check\n");
    fprintf(stderr,"\t-d      decompress\n");
    fprintf(stderr,"\t-i ext  extension for the input file  (def. .bin)\n");
    fprintf(stderr,"\t-o ext  extension for the output file (def. .k2x)\n");
    fprintf(stderr,"\t-m M    minimatrix size (def. 2), compression only\n");
    fprintf(stderr,"\t-v      verbose\n\n");
    fprintf(stderr,"Default action is to compress basename.bin to basename.k2x\n"
    "To decompress use -d option and appropriate extensions eg \"-d -i .k2x -o .k2d\"\n\n");
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

