/* Multiplication of boolean matrices represented in bbm format
   using opm 
    
   In medium size experiments it was sligthly faster than the serial version
   1 secs vs 2 for a 18k matrix.  
    
   Note that sparsity does affect the running time

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
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>
#include "bbm.h"
#include <omp.h>

#define default_ext ".prod"

static void usage_and_exit(char *name);
static void quit(const char *msg, int line, char *file);

static void fast_mmult_bbm(const uint8_t *a, size_t size, const uint8_t *b, uint8_t *c);


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
  int threads = omp_get_max_threads();
  char *ext = default_ext;
  while ((c=getopt(argc, argv, "e:cnvt:")) != -1) {
    switch (c) 
      {
      case 'e':
        ext = optarg; break;                 
      case 't':
        threads = atoi(optarg); break;                 
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
    fputs("==== Command line:\n",stdout);
    for(int i=0;i<argc;i++)
     fprintf(stdout," %s",argv[i]);
    fputs("\n",stdout);  
  }

  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 3) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  fprintf(stdout,"Working with %d threads out of %d\n",
          threads,omp_get_max_threads());
  omp_set_num_threads(threads);


  // create file names
  sprintf(iname1,"%s",argv[1]);
  sprintf(iname2,"%s",argv[2]);
  sprintf( oname,"%s%s",argv[1],ext); 

  // read matrices
  uint8_t *a, *b, *ab;
  size_t size;
  a = bbm_read(iname1,&size);
  if(strcmp(iname1,iname2)==0)
    b = a;
  else {
    size_t size1;
    b = bbm_read(iname2,&size1);
    if(size1!=size) quit("Input matrices have different sizes",__LINE__,__FILE__);
  }
  // compute product
  ab = bbm_alloc(size);
  time_t start_mult = time(NULL);
  fast_mmult_bbm(a,size,b,ab);
  fprintf(stderr,"fast_mmult_bbm elapsed time: %.0lf secs\n",(double) (time(NULL)-start_mult));
  if(write) bbm_write(ab,size,oname);

  if(check) {
    uint8_t *m4 = bbm_alloc(size);
    start_mult = time(NULL);
    mmult_bbm(a,size,b,m4);
    fprintf(stderr,"mmult_bbm elapsed time: %.0lf secs\n",(double) (time(NULL)-start_mult));
    ssize_t eq = mequals_bbm(ab,size,m4);
    if(eq<0) fprintf(stdout,"Product matches the one computed using standard algorithm!\n");
    else fprintf(stdout,"Product matrix differs at position (%zd,%zd) "
      "fast:%d vs std:%d\n",eq/size,eq%size, ab[eq], m4[eq]);
    free(m4);    
  }

  // free and terminate
  free(a);
  if(strcmp(iname1,iname2)) free(b);
  free(ab);    
  // report running time
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));
  if(verbose) fprintf(stderr,"==== Done\n");
  return EXIT_SUCCESS;
}


void fast_mmult_bbm(const uint8_t *a, size_t size, const uint8_t *b, uint8_t *c) {
  assert(a!=NULL && b!=NULL && c!=NULL && size>0);
  // clean c
  byte_to_bbm(c,size,0,0,size,0);

  #pragma omp parallel for          // this is currently the only parallelization
  for(size_t i=0; i<size; i++)
    for(size_t k=0; k<size; k++) 
      if(a[i*size+k]) {
        for(size_t j=0; j<size; j++) { // note: this for also is parallelizable
          if(b[k*size+j])
            c[i*size+j] = 1;
        }
      }
}




static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] iname1 iname2\n\n",name);
    fputs("Options:\n",stderr);
    fprintf(stderr,"\t-t thr  number of omp threads (def. %d)\n",omp_get_max_threads());
    fprintf(stderr,"\t-n      do not write output file, only show stats\n");    
    fprintf(stderr,"\t-e ext  extension for the output file (def. %s)\n",default_ext);
    fprintf(stderr,"\t-c      check multiplication (very slow for large matrices!)\n");
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

