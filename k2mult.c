/* Multiplication of boolean matrices represented as k2 trees

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
  char iname1[PATH_MAX], iname2[PATH_MAX], oname[PATH_MAX];
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  bool check = false;
  char *iext = "", *oext=NULL;
  while ((c=getopt(argc, argv, "i:o:cv")) != -1) {
    switch (c) 
      {
      case 'i':
        iext = optarg; break;  
      case 'o':
        oext = optarg; break;                
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
  if (oext==NULL) oext = iext; // default oext is iext


  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 4) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // create file names
  sprintf(iname1,"%s%s",argv[1],iext);
  sprintf(iname2,"%s%s",argv[2],iext);
  sprintf( oname,"%s%s",argv[3],oext); 

  // init k2 library
  k2mat_t a = K2MAT_INITIALIZER, b=K2MAT_INITIALIZER, ab=K2MAT_INITIALIZER;
  int size, asize;

  size = mload(&asize, &a, iname1);
  if (verbose) show_stats(asize,&a,"A");
  if(strcmp(iname1,iname2)==0) {
    b = a;
    b.read_only = true;
  } else {
    int bsize, size1 = mload(&bsize, &b, iname2);
    if(size1!=size) quit("Input matrices have different sizes",__LINE__,__FILE__);
    if(bsize!=asize) quit("k2 matrices have different sizes",__LINE__,__FILE__);
  }
  if (verbose) show_stats(asize,&b,"B");
  mmult(asize,&a,&b,&ab);
  if (verbose) show_stats(asize,&ab,"A*B");
  msave(size,asize,&ab,oname);
  k2_free(&a);
  if(strcmp(iname1,iname2)) k2_free(&b);
  k2_free(&ab);    
  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  (void) check;
  return EXIT_SUCCESS;
}


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] iname1 iname2 oname\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-c      check multiplication\n");
    fprintf(stderr,"\t-d      decompress\n");
    fprintf(stderr,"\t-i ext  extension for the input files  (def. none)\n");
    fprintf(stderr,"\t-o ext  extension for the output file (def. same as input)\n");
    fprintf(stderr,"\t-v      verbose\n\n");
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

