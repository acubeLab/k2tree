/* Testing of arithmetic/logical operations on binary matrices 
   represented as k^2 trees 

   Test of the operations in k2ops.c

   Note: internally matrix dimensions are always of the form 2^k times the size 
   of a minimatrix (those stored at the leaves of the tree), with k>0, 
   but the input can be of any size, and the matrix will be padded with 0's 

   Technical note: matrix sizes are int, therefore limited to 2^31, but values related to
   the overall number of elements is always stored into a size_t variable (usually 64 bits))  

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

// extension for boolean byte matrix format
#define File_ext ".bbm"

static void usage_and_exit(char *name);
// static void quit(const char *msg, int line, char *file);
// static int intsqrt(int n);



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
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  int mmsize = 2;
  while ((c=getopt(argc, argv, "m:v")) != -1) {
    switch (c) 
      {
      case 'm':
        mmsize = atoi(optarg); break;
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

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;


  // init k2 library
  minimat_init(mmsize);
  // read matrix from file
  size_t size;
  uint8_t *buffer = bbm_read(argv[1],&size);
  if(verbose>0) bbm_to_ascii(buffer,size,0,0,size,stderr);
  // convert to k2mat_t
  k2mat_t a = K2MAT_INITIALIZER, b = K2MAT_INITIALIZER;
  int asize = mread_from_bbm(buffer,size,&a);
  mwrite_to_bbm(buffer,size,asize,&a);
  if(verbose>0) bbm_to_ascii(buffer,size,0,0,asize,stderr);
  show_stats(asize,&a,"A");
  mmult(asize ,&a,&a,&b); // b = a*a
  mwrite_to_bbm(buffer,size,asize,&b);
  if(verbose>0) bbm_to_ascii(buffer,size,0,0,asize,stderr);
  show_stats(asize,&b,"A*A");
   
  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  
  return EXIT_SUCCESS;
}

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] infile\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-m M    minimatrix size (def. 2)\n");
    fprintf(stderr,"\t-v      verbose\n\n");
    exit(1);
}

/*
// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

// compute integer square root
static int intsqrt(int n) {
  assert(n>=0);
  int x = n;
  int y = (x + 1) / 2;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2;
  }
  assert(x*x <= n && (x+1)*(x+1) > n);
  return x;
} 
*/
