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

#define default_ext "prod.k2"

static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);


int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  char iname1[PATH_MAX], iname2[PATH_MAX], oname[PATH_MAX];
  time_t start_wc = time(NULL);

  /* ------------- read options from command line ----------- */
  opterr = 0;
  bool check = false, write = true;
  char *ext = default_ext;
  while ((c=getopt(argc, argv, "e:cnv")) != -1) {
    switch (c) 
      {
      case 'e':
        ext = optarg; break;                 
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


  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 3) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // create file names
  sprintf(iname1,"%s",argv[1]);
  sprintf(iname2,"%s",argv[2]);
  sprintf( oname,"%s%s",argv[1],ext); 

  // init k2 variables
  k2mat_t a = K2MAT_INITIALIZER, b=K2MAT_INITIALIZER, ab=K2MAT_INITIALIZER;
  int size, asize;

  size = mload_from_file(&asize, &a, iname1); // also init k2_library
  if (verbose) mshow_stats(asize,&a,iname1,stderr);
  if(strcmp(iname1,iname2)==0) {
    k2make_pointer(&a,&b);
  } else {
    int bsize, size1 = mload_from_file(&bsize, &b, iname2);
    if(size1!=size) quit("Input matrices have different sizes",__LINE__,__FILE__);
    if(bsize!=asize) quit("k2 matrices have different sizes",__LINE__,__FILE__);
  }
  if (verbose) mshow_stats(asize,&b,iname2,stderr);
  mmult(asize,&a,&b,&ab);
  if (verbose || !write) mshow_stats(asize,&ab,oname,stderr);
  if(write) msave_to_file(size,asize,&ab,oname);

  // check product if requested
  if(check) {
    uint8_t *m1 = bbm_alloc(size), *m2 = bbm_alloc(size), *m3 = bbm_alloc(size); 
    mwrite_to_bbm(m1,size,asize,&a);
    mwrite_to_bbm(m2,size,asize,&b);
    mwrite_to_bbm(m3,size,asize,&ab);
    uint8_t *m4 = bbm_alloc(size);
    mmult_bbm(m1,size,m2,m4);
    ssize_t eq = mequals_bbm(m3,size,m4);
    if(eq<0) fprintf(stderr,"Product is correct!\n");
    else fprintf(stderr,"Product matrix differs at position (%zd,%zd) "
      "k2:%d vs bbm:%d\n",eq/size,eq%size, m3[eq], m4[eq]);
    free(m1); free(m2); free(m3); free(m4);
  }

  // free and terminate
  k2_free(&a);
  if(strcmp(iname1,iname2)) k2_free(&b);
  k2_free(&ab);    
  // report running time
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  fprintf(stderr,"==== Done\n");
  (void) check;
  return EXIT_SUCCESS;
}


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] iname1 iname2\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-c      check multiplication\n");
    fprintf(stderr,"\t-n      do not write output file, only show stats\n");    
    fprintf(stderr,"\t-e ext  extension for the output file (def. %s)\n",default_ext);
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

