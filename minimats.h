/* Routines for arithmetic/logical operations on binary matrices 
   represented as k^2 trees 

   This file contains the fucntion releated to minimatrices, ie,
   the matrices stored as leaves of the k^2 tree.

   minimat matrices have size MMsize which is set (and never changed)
   at the start of the execution 

   Copyright August 2023-today   ---  giovanni.manzini@unipi.it
*/

static void quit(const char *msg, int line, char *file);

// the type representing a (temporary) minimat matrix must be an integer
// so that two minimat matrices can be summed with a single
// bitwise OR operation
typedef uint64_t minimat_t; // explict matrix aka minimat (currently 2x2)

// ----------------------------------------------
// global variable storing the size of a mini-matrix 
// ie. the last level of recursion
// possibly this will be a command line parameter
static int MMsize = 2; 
// global variable storing how many nodes takes a minimat:
// since nodes have 4 bits (not likely to change)
// this amounts to how many nibbles takes a minimat 
// a 2x2 binary matrix takes one nibble, so here the constant is 1
static const int Minimat_node_ratio = 1;
#if ((MMsize*MMsize) != 4*Minimat_node_ratio)
#error "MMsize vs Minimat_node_ratio mismatch"
#endif

// minimat constants (depend on size, here are 2x2)
#define MINIMAT0s      0x0   // minimat containing all 0's  
#define MINIMAT1s      0xF   // minimat containing all 1's  

// check that minimat_t is large enough to store a minimat
#if ((MMsize*MMsize) > 64 ) // 64=8*sizeof(minimat_t), update if minimat_t changes
#error "MMsize too large for minimat_t"
#endif



// global variable containing the product of every pair of possible minimatrices
// to initialized in minimats_init();
static minimat_t *mprods2x2 = NULL; // will contain 256 entries  ;

// multiply two minimat matrices of size 2*2
minimat_t mmult2x2(minimat_t a, minimat_t b) {
  return mprods2x2[(a<<4) | b]; // equivalent to minimat_prods[a][b]
}

void init_mprods2x2(void) {
  // create array
  assert(mprods2x2==NULL);
  mprods2x2 = malloc(256*sizeof(minimat_t));
  if(mprods2x2==NULL) quit("init_mprods2x2: malloc failed",__LINE__,__FILE__);
  // fill with products
  minimat_t arow[2], bcol[2],c;
  for(minimat_t a=0; a<16; a++) 
    for(minimat_t b=0; b<16;b++) {
      arow[0] = a&3; arow[1] = a>>2;
      bcol[0] = (b&1) | (b&4)>>1; 
      bcol[1] = (b&2)>>1 | (b&8)>>2;
      c =  (arow[0] | bcol[0]) ? 1 : 0; // c[0][0] = a[0]*b[0]
      c |= (arow[0] | bcol[1]) ? 2 : 0; // c[0][1] = a[0]*b[1]
      c |= (arow[1] | bcol[0]) ? 4 : 0; // c[1][0] = a[1]*b[0]
      c |= (arow[1] | bcol[1]) ? 8 : 0; // c[1][1] = a[1]*b[1]
      mprods2x2[a<<16 | b] = c;
    }
}

void minimats_init(int msize) {
  if(MMsize!=2) quit("minimats_init: MMsize!=2",__LINE__,__FILE__);
  MMsize = msize;
  init_mprods2x2();
}


// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);

  exit(1);
}
