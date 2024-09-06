/* Pagerank for graph matrices represented as k2 trees

   Note: internally matrix dimensions are always of the form 2^k times the size 
   of a minimatrix (those stored at the leaves of the tree), with k>0, 
   but the input can be of any size, and the matrix will be padded with 0's 
  
 * Overview
 * In the original formulation PageRank requires the left multiplication 
 * of the current rank vector, with entries divided by the outdegree
 * of each node, times the binary adjacency matrix.
 * 
 * Since we have implemented only right matrix-vector multiplication,
 * we are assuming that the graph matrix has been 
 * transposed (and all self-loops already removed).  
 * 
 * The main iteration goes as follows:
 *   1 each rank_i entry is divided by the outdegree of node i
 *     which is now the # of nonzero elements in column i
 *   2 if column i has no nonzero then i is a dandling node 
 *     and rank_i is added to the sum of dandling nodes ranks
 *   3 the normalized rank vector is right multiplied by the 
 *     (transpose adjacency matrix)
 *   4 the new rank vector is obtained for the above product
 *     + contribution of teleporting and dandling nodes.
 * In mathematical terms, let
 *   N = # nodes
 *   d = damping factor (from the command line)
 *   X = current rank vector
 *   dnr = 0  # dandling nodes rank sum
 *   for i in range(N):
 *     if col_count[i]==0: dnr[i] += X[i]
 *     else Y[i] = X[i]/col_count[i]
 *   Z = M*Y
 *   for i in range(N):
 *     Z[i] = d*Z[i] + (d/N)*dnr + (1-d)/N
 *   X = Z  # prepare for next iteration    
 * 
 * We start the iterations with X = (1/N ... 1/N)^T
 * we stop after maxiter (from the command line) iterations or
 * when the sum of the abs differences of the ranks between two
 * consecutive iterations is smaller than eps (from the command line)   
 *  
 * The above is the standard algorithm; in this program we use only two 
 * vectors X and Z at the expense of some loss of accuracy in the 
 * error computation. The trick is that in X and Y there is the same 
 * information provided we know the outdegree of each node; so
 * we do not store X but only Y and retrieve X from Y when needed
 *   1. compute Y overwriting X and if col_count[i]==0 simply set Y[i] = X[i]
 *   2. compute Z = M*Y
 *   3. compute the new rank values but do not store them, instead
 *      compute the error, dnr, and Y for the next iteration
 
 * Note: the program can be compiled with -DUSE_BARRIER to use pthread barriers
 *       instead of semaphores for synchronization between the main thread and
 *       the worker threads. Experiments show that the barrier version is
 *       minimally faster than the semaphore version (183 vs 185 for it-2004) 
 *       using 8 threads on a machine with plenty of cores. 
 *       Compiling with -DDETAILED_TIMING provides the time spent in matrix
 *       multiplication only which turns out to use roughly 90% of the total time
 *       (168 out of 184 for it-2004 8 threads) this suggests that 
 *       parallelizing the vector update and error computation beween 
 *       multiplications would not provide a significant speedup.
 *       Similar results are obtained for arabic-2005.
 *       Experiments for 16 threads on it-2004 confirm this trend:
 *          semaphores:  elapsed 107, multiplication only: 89 secs
 *          barriers:    elapsed 106, multiplication only: 90 secs
 *       I tried to measure the memory usage with malloc_count but it seems
 *       that the library is no longer working properly.
 * 
 *       Although the code contains macros referring to the b128mat econding
 *       it cannot be compiled with -DB128MAT because the b128 library is
 *       missing (at least) the matrix vector multiplication functions.   
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
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <sys/times.h>
// definitions to be used for b128 vs k2-encoded matrices 
#ifdef B128MAT
#include "b128.h"
#define K2MAT_INITIALIZER B128MAT_INITIALIZER
typedef b128mat_t k2mat_t;
#else // k2mat
#include "k2.h"
#endif
// used by both matrix type
#include "bbm.h"
#include "extra/xerrors.h"



// for compatibilty with matrepair
typedef struct {
  size_t size; // size of the vector
  vfloat *v; // vector
} vector;


static void vector_destroy(vector *v)
{
  if(v->v!=NULL) free(v->v);
  free(v);
}

// input/output data for each thread 
typedef struct {
  k2mat_t *a;    // matrix block (but with same size whole martix)
  vector *rv;    // right vector  
  vector *lv;    // left vector
  size_t size;   // size of the input/output vectors
  size_t asize;  // internal size of the matrix block
  int op;        // operation to execute
  #ifndef USE_BARRIER
  sem_t *in;     // semaphore for input shared with the main thread
  sem_t *out;    // semaphore for output shared with the main thread
  #else
  pthread_barrier_t *barrier; // barrier for synchronization
  #endif
} tdata;




// static functions at the end of the file
static void quit(const char *msg, int line, const char *file);
static void minHeapify(double v[], int arr[], int n, int i);
static void kLargest(double v[], int arr[], int n, int k);
static vector *vector_create_value(size_t n, vfloat v);
static void vector_destroy(vector *v);
static void mload_from_file_multipart(size_t *, size_t, k2mat_t *, int, char *, char *);
static void *block_main(void *v);


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] matrix col_count_file\n",name);
    fprintf(stderr,"\t\t-v             verbose\n");
    fprintf(stderr,"\t\t-b num         number of row blocks, default: 1\n");
    fprintf(stderr,"\t\t-x ext         extension (only added when b>1), default: \".k2\"\n");
    fprintf(stderr,"\t\t-m maxiter     maximum number of iterations, default: 100\n");
    fprintf(stderr,"\t\t-e eps         stop if error<eps (default: ignore error)\n");
    fprintf(stderr,"\t\t-d df          damping factor (default: 0.9)\n");
    fprintf(stderr,"\t\t-k K           show top K nodes (default: 3)\n\n");
    exit(1);
}

int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int c;
  time_t start_wc = time(NULL);
  #ifdef DETAILED_TIMING
  struct tms ignored;
  clock_t t1,t2;
  long m1=0;
  #endif
  // default values for command line parameters 
  int maxiter=100,nblocks=1,topk=3;
  double dampf = 0.9, eps = -1;
  char *ext = ".k2";
  
  /* ------------- read options from command line ----------- */
  opterr = 0;
  while ((c=getopt(argc, argv, "m:e:d:k:vb:x:")) != -1) {
    switch (c) 
      {
      case 'v':
        verbose++; break;
      case 'm':
        maxiter=atoi(optarg); break;
      case 'e':
        eps=atof(optarg); break;
      case 'd':
        dampf=atof(optarg); break;
      case 'k':
        topk=atoi(optarg); break;
      case 'b':
        nblocks=atoi(optarg); break; 
      case 'x':
        ext=optarg; break;
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
  // check command line
  if(maxiter<1 || topk<1 || nblocks<1) {
    fprintf(stderr,"Error! Options -b -m and -k must be at least one\n");
    usage_and_exit(argv[0]);
  }
  if(dampf<0 || dampf>1) {
    fprintf(stderr,"Error! Options -d must be in the range [0,1]\n");
    usage_and_exit(argv[0]);
  }
  
  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 3) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;

  // ----------- read column count file and get matrix size 
  u_int32_t *outd = NULL;
  size_t size;
  {
    FILE *ccol_file  = fopen(argv[2],"rb");
    if(ccol_file==NULL) quit("Cannot open col_count_file", __LINE__, __FILE__);
    long e = fseek(ccol_file,0,SEEK_END);
    if(e!=0) quit("fseek failed on col_count_file", __LINE__, __FILE__);
    e = ftell(ccol_file)/sizeof(u_int32_t);
    if(e<1) quit("ftell failed or invalid col_count_file", __LINE__, __FILE__);
    size = (size_t) e;
    outd = (u_int32_t *) malloc(size*sizeof(*outd));
    if(outd==NULL) quit("Cannot allocate out_degree vector", __LINE__, __FILE__);
    rewind(ccol_file);
    e = fread(outd,sizeof(*outd),size,ccol_file);
    if(e!=size) quit("cannot read out_degree vector from col_count_file", __LINE__, __FILE__);
    if(fclose(ccol_file)!=0) quit("Error closing col_count_file", __LINE__, __FILE__);
  }

  // compute # dandling nodes and arcs if verbose
  if(verbose>0) {
    fprintf(stderr,"Number of nodes: %ld\n",size);
    long dn=0,arcs=0;
    for(int i=0;i<size;i++)
      if(outd[i]==0) dn++;
      else arcs += outd[i];
    fprintf(stderr,"Number of dandling nodes: %ld\n",dn);
    fprintf(stderr,"Number of arcs: %ld\n",arcs);
  }

  // init k2 variables: single matrix (nblock=1) 
  k2mat_t a=K2MAT_INITIALIZER;
  size_t asize;
  // init k2 matrices nblocks>1
  // blocks are sets of rows so that matrix multiplication can be parallelized
  // each block is acutally a full matrix with the same size of the original matrix
  k2mat_t rblocks[nblocks];  
  for(int i=0;i<nblocks;i++) rblocks[i]=a; // struct with all fields set to 0 or NULL

  if(nblocks==1) {
    size_t msize = mload_from_file(&asize, &a, argv[1]); // also init k2 library
    if(msize!=size) quit("Matrix size mismatch", __LINE__, __FILE__);
    if (verbose>1) mshow_stats(size,asize,&a,argv[1],stdout);
  }
  else {
    mload_from_file_multipart(&asize, size, rblocks, nblocks, argv[1], ext);
    if (verbose>1) {
      for(int i=0;i<nblocks;i++) {
        char name[FILENAME_MAX];
        sprintf(name,"%s.%d.%d%s",argv[1],nblocks,i,ext);
        mshow_stats(size,asize,&rblocks[i],name,stdout);
      }
    }
  }

  // ------------ alloc rank and aux vectors
  vector *y = vector_create_value(size,0);
  vector *z = vector_create_value(size,0);

  // data structures for multithread computation (nblocks>1)
  tdata td[nblocks];
  pthread_t t[nblocks];
  #ifndef USE_BARRIER
  sem_t tsem_in, tsem_out;
  #else
  pthread_barrier_t tbarrier;
  #endif

  // initialize thread data and start threads
  if(nblocks>1) {
    #ifndef USE_BARRIER
    xsem_init(&tsem_in,0,0,__LINE__,__FILE__);
    xsem_init(&tsem_out,0,0,__LINE__,__FILE__);
    #else
    pthread_barrier_init(&tbarrier,NULL,nblocks+1);
    #endif
    for(int i=0;i<nblocks;i++) {
      td[i].a = &rblocks[i];
      td[i].size = size;
      td[i].asize = asize;
      td[i].rv = y;      // right vector
      td[i].lv = z;      // left vector
      td[i].op = 0;      // right multiplication
      #ifndef USE_BARRIER
      td[i].in = &tsem_in;
      td[i].out = &tsem_out;
      #else
      td[i].barrier = &tbarrier;
      #endif
      xpthread_create(&t[i],NULL,&block_main,&td[i],__LINE__,__FILE__);
    }
  }

  // x_0 = (1/N ... 1/N)^T (no need to write those values)
  // initialze y_0 based on x_0=(1/N ... 1/N)^T and compute dandling nodes rank
  double dnr = 0;
  for(int i=0;i<size;i++) 
    if(outd[i]==0) dnr += (y->v[i] = 1.0/(double)size);
    else y->v[i] = (1.0/(double)size)/outd[i];
  // main loop 
  int iter=0;
  double delta=11+eps; // always larger than eps to prevent stopping at first iteration
  while(iter<maxiter && delta>=eps) {
    #ifdef DETAILED_TIMING
    t1 = times(&ignored);
    #endif
    if (nblocks==1) 
      mvmult(asize,&a,size,y->v,z->v, true);       // z = M*y
    else {
      #ifndef USE_BARRIER
      for(int i=0;i<nblocks;i++) xsem_post(&tsem_in,__LINE__,__FILE__);   // start threads
      for(int i=0;i<nblocks;i++) xsem_wait(&tsem_out,__LINE__,__FILE__);  // wait threads 
      #else
      pthread_barrier_wait(&tbarrier);
      pthread_barrier_wait(&tbarrier);
      #endif
    }
    #ifdef DETAILED_TIMING
    t2 = times(&ignored);
    m1 += (t2-t1);       // measure time for matrix multiplication only
    #endif
    // compute contribution of teleporting and dandling nodes
    double teleport = (dnr*dampf+1-dampf)/(double)size;
    // compute new rank vector values (without storing them) and
    //   1. the L1 norm of the difference with the previos iteration retrifred from y
    //   2. update y with the new values, and compute the new dnr
    dnr = delta = 0;
    for(int i=0;i<size;i++) {
      double nextri = dampf*z->v[i] + teleport;
      z->v[i]=0; // clear z for next multiplication operation
      if (outd[i]==0) {
        delta += fabs(nextri-y->v[i]); // update delta
        dnr += (y->v[i]=nextri);       // update dnr and y_i
      }
      else {
        delta += fabs(nextri-y->v[i]*outd[i]); //update delta
        y->v[i] = nextri/outd[i];              // update y_i
      }
    }
    iter++; // iteration complete
    if(verbose>1) fprintf(stderr,"Iteration %d, delta=%g \n",iter,delta);
  }
  // stop and join threads
  if(nblocks>1) {
    for(int i=0;i<nblocks;i++) {
      td[i].op = -1;
      #ifndef USE_BARRIER 
      xsem_post(&tsem_in,__LINE__,__FILE__);
      #endif
    }
    #ifdef USE_BARRIER
    pthread_barrier_wait(&tbarrier);
    #endif
    for(int i=0;i<nblocks;i++)
      xpthread_join(t[i],NULL,__LINE__,__FILE__);
  }
  // retrieve the actual rank vector from the last iteration
  for(int i=0;i<size;i++) 
    if(outd[i]!=0) y->v[i] = y->v[i]*outd[i];
  // call x the actual rank vector for consistency with the literature and disable y
  vector *x = y; y = NULL; // make sure y is never used again

  if(verbose>0) {
    if (delta>eps) fprintf(stderr,"Stopped after %d iterations, delta=%g\n",iter,delta);
    else           fprintf(stderr,"Converged after %d iterations, delta=%g\n",iter,delta);
    double sum=0;
    for(int i=0;i<size;i++) sum += x->v[i];
    fprintf(stderr,"Sum of ranks: %f (should be 1)\n",sum);
  }  
  // free k2 matrices
  if(nblocks==1) matrix_free(&a);
  else for(int i=0;i<nblocks;i++) matrix_free(&rblocks[i]);
  minimat_reset(); // reset the minimat library and free minimat product table
  // deallocate z and outd: we may need space for the topk array
  vector_destroy(z);
  free(outd);

  // retrieve topk nodes
  if(topk>size) topk = (int) size;
  int *top = (int *) malloc(topk*sizeof(*top));
  int *aux = (int *) malloc(topk*sizeof(*top));
  if(top==NULL || aux==NULL) quit("Cannot allocate topk/aux array", __LINE__, __FILE__);
  kLargest(x->v,aux,(int)size,topk); // only work for size <2^31, fix this!!!
  // get sorted nodes in top
  for(int i=topk-1;i>=0;i--) {
    top[i] = aux[0];
    aux[0] = aux[i];
    minHeapify(x->v,aux,i,0);
  }  
  // report topk nodes sorted by decreasing rank
  if (verbose>0) {
    fprintf(stderr, "Top %d ranks:\n",topk);
    for(int i=0;i<topk;i++) fprintf(stderr,"  %d %lf\n",top[i],x->v[top[i]]);
  }
  // report topk nodes id's only on stdout
  fprintf(stdout,"Top:");
  for(int i=0;i<topk;i++) fprintf(stdout," %d",top[i]);
  fprintf(stdout,"\n");
  // destroy everything
  free(top); free(aux);
  vector_destroy(x);
  #ifdef DETAILED_TIMING
  fprintf(stderr,"Total mult time (secs): %lf  Average: %lf\n", ((double)m1)/((double)sysconf(_SC_CLK_TCK)), 
                                                         ((double)m1/(double)iter)/((double) sysconf(_SC_CLK_TCK)));
  #endif
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));  
  return 0;
}


// function executed by each thread 
// wait on a semaphore for a new operation to execute
// on its given matrix
static void *block_main(void *v)
{
  tdata *td = (tdata *) v;
  assert(td->rv!=NULL); // 
  assert(td->lv!=NULL); //  
  assert(td->size <= td->asize);
  while(true) {
    // wait for input
    #ifndef USE_BARRIER 
    xsem_wait(td->in,__LINE__,__FILE__);
    #else
    pthread_barrier_wait(td->barrier);
    #endif
    if(td->op<0) break;  // exit loop
    else if(td->op==0) { //right mult
      mvmult(td->asize,td->a,td->size,td->rv->v,td->lv->v,false);
    }
    else quit("Unknown operation", __LINE__, __FILE__);
    // output ready
    #ifndef USE_BARRIER 
    xsem_post(td->out,__LINE__,__FILE__);
    #else
    pthread_barrier_wait(td->barrier);
    #endif
  }
  return NULL;
}


// load the nb blocks of a matrix and store them in the array of k2mat_t
static void mload_from_file_multipart(size_t *asize, size_t size, k2mat_t *b, int nb, char *name, char *ext) {
  assert(nb>1);
  size_t as=0, msize; 
  for(int i=0;i<nb;i++) {
    char fname[FILENAME_MAX];
    sprintf(fname,"%s.%d.%d%s",name,nb,i,ext);
    msize = mload_from_file(asize, &b[i], fname); // also init k2 library
    if(msize!=size) quit("Matrix size mismatch", __LINE__, __FILE__);
    if(as==0) as = *asize;
    else if(as!=*asize) quit("Internal matrix size mismatch", __LINE__, __FILE__);
  }
}


// heap based algorithm for finding the k largest ranks in heap order

// A utility function to swap two elements
static void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

// Heapify a subtree rooted with node i which is an index in arr[]. 
// n is size of heap. the key associated to entry arr[i] is v[arr[i]]  
static void minHeapify(double v[], int arr[], int n, int i) {
    int smallest = i;  // Initialize smallest as root
    int left = 2*i + 1;
    int right = 2*i + 2;

    // If left child is smaller than root
    if (left < n && v[arr[left]] < v[arr[smallest]])
        smallest = left;

    // If right child is smaller than smallest so far
    if (right < n && v[arr[right]] < v[arr[smallest]])
        smallest = right;

    // If smallest is not root
    if (smallest != i) {
        swap(&arr[i], &arr[smallest]);
        // Recursively heapify the affected sub-tree
        minHeapify(v, arr, n, smallest);
    }
}

// Function to find the k'th largest elements in an array
// v[0..n-1], arr[0..k-1] is the output array already allocated
static void kLargest(double v[], int arr[], int n, int k) {
  assert(k<=n);
  assert(k>0);
  // init arr[] with the first k elements
  for(int i=0;i<k;i++) arr[i] = i;
  // Build a min heap of first (k) elements in arr[]
  for (int i = k / 2 - 1; i >= 0; i--)
    minHeapify(v, arr, k, i);
  // Iterate through the rest of the array elements
  for (int i = k; i < n; i++) {
    // If current element is larger than root of the heap
    if (v[i] > v[arr[0]]) {
      // Replace root with current element
      arr[0] = i;
      // Heapify the root
      minHeapify(v, arr, k, 0);
    }
  }
}



static vector *vector_create_value(size_t n, vfloat v)
{
  vector *w = malloc(sizeof(vector));
  if(w==NULL) quit("malloc failed", __LINE__, __FILE__);
  w->size=n; w->v = (vfloat *) malloc(n*sizeof(vfloat));
  if(w->v==NULL) quit("malloc failed", __LINE__, __FILE__);
  for(size_t i=0;i<n;i++) w->v[i]=v;
  return w;
}

// write error message and exit
static void quit(const char *msg, int line, const char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

