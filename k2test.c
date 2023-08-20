/* Testing of arithmetic/logical operations on binary matrices 
   represented as k^2 trees 

   Testo of the operations in k2ops.c

   Matrix dimensions are assumed to be power of 2, of size at least
   2*MMSize (minimatrix size), ie, the size of the last level of recursion. 


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
#include "k2.h"


// extension for boolean byte matrix format
#define File_ext ".bbm"



static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] leftmatrix rightmatrix outmatrix\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-v             verbose\n\n");
    exit(1);
}

// write error message and exit
void quit(const char *s)
{
  if(errno==0) fprintf(stderr,"%s\n",s);
  else  perror(s);
  exit(1);
}

// compute integer square root
int intsqrt(int n) {
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



int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  while ((c=getopt(argc, argv, "b:cfinv")) != -1) {
    switch (c) 
      {
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
  if (argc-optind != 1) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  minimats_init(2); // minimatrix size = 2x2
  k2mat_t a = K2MAT_INITIALIZER, b = K2MAT_INITIALIZER;
  int size =4;
  mmult(size ,&a,&a,&b); // b = a*a
  size_t pos, nodes, minimats;
  int levels = mstats(size,&b,&pos,&nodes,&minimats);
  printf("Levels: %d, Pos: %zd, Nodes: %zd, Minimats: %zd\n",levels,pos,nodes,minimats);
  
  // ----------- read and check # rows and cols 
  // statistics
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  
  return EXIT_SUCCESS;
}
