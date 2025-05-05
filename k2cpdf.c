// c includes
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <string.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <stdbool.h>
#include <libgen.h>
#include <inttypes.h>
#include <stdio.h> 

// local includes
#include "k2.h"
#include "libsais/include/libsais64.h"

// from succint repository
const uint8_t debruijn64_mapping[64] = {
  63,  0, 58,  1, 59, 47, 53,  2,
  60, 39, 48, 27, 54, 33, 42,  3,
  61, 51, 37, 40, 49, 18, 28, 20,
  55, 30, 34, 11, 43, 14, 22,  4,
  62, 57, 46, 52, 38, 26, 32, 41,
  50, 36, 17, 19, 29, 10, 13, 21,
  56, 45, 25, 31, 35, 16,  9, 12,
  44, 24, 15,  8, 23,  7,  6,  5
};

const uint64_t debruijn64 = 0x07EDD5E59A4E28C2ULL;

uint8_t bit_position(uint64_t x){
    return debruijn64_mapping[(x * debruijn64) >> 58];
}

uint8_t _msb(uint64_t x, unsigned long* ret){
  if (!x)
    return false;

  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;

  x ^= x >> 1;
  *ret = bit_position(x);

  return true;
}

uint8_t msb(uint64_t x){
  unsigned long ret = -1U;
  _msb(x, &ret);
  return (uint8_t)ret;
}

uint64_t ceil_log2(uint64_t x) {
  return (x > 1) ? msb(x - 1) + 1 : 0;
}

uint64_t floor_log2(uint64_t x) {
  return (x > 1) ? msb(x) : 0;
}

typedef struct {
  size_t k, maxn;
  uint64_t **st;
} st_t;

int64_t min(int64_t a, int64_t b) {
  return a < b ? a : b;
}

void st_init(st_t *st, size_t n, int64_t *arr) {
  st->k = floor_log2(n) + 1;
  st->maxn = n + 1;
  st->st = malloc(sizeof(uint64_t*) * (st->k + 1));
  for(size_t i = 0; i < st->k; i++) {
    st->st[i] = calloc(st->maxn, sizeof(uint64_t));
  }

  for(size_t i = 0; i < st->maxn; i++) {
    st->st[0][i] = arr[i];
  }

  for(size_t j = 1; j < st->k; j++) {
    for(size_t i = 0; i + (1 << j) <= n; i++) {
      st->st[j][i] = min(st->st[j - 1][i], st->st[j - 1][i + (1 << (j - 1))]);
    }
  }
}

uint64_t st_query(st_t *st, size_t l, size_t r) {
  size_t k = floor_log2(r - l + 1);
  return min(st->st[k][l], st->st[k][r - (1 << k) + 1]);
}

void st_free(st_t *st) {
  for(size_t i = 0; i < st->k; i++) {
    free(st->st[i]);
  }
  free(st->st);
}

typedef struct {
  uint64_t n;
  int64_t* e;
} dsu;

void dsu_init(dsu* u, uint64_t n) {
  u->e = (int64_t*) malloc(sizeof(int64_t) * n);
  u->n = n;
  for(size_t i = 0; i < n; i++) u->e[i] = -1;
}

int64_t dsu_find_set(dsu* u, int64_t x) {
  if(u->e[x] < 0) return x;
  u->e[x] = dsu_find_set(u, u->e[x]);
  return u->e[x];
}

int8_t dsu_union_set(dsu* u , int64_t x, int64_t y) {
  x = dsu_find_set(u, x);
  y = dsu_find_set(u, y);
  if(x == y) return 0;
  if(x > y) {
    int64_t aux = x;
    x = y;
    y = aux;
  }
  u->e[x] += u->e[y];
  u->e[y] = x;
  return 1;
}

void dsu_free(dsu* u) {
  free(u->e);
}

uint64_t get_size_rec(uint8_t *tree, vu64_t *z, uint64_t t_pos, uint64_t t_size,
                      uint64_t z_pos, uint64_t b_tree) {
  uint64_t c_root = __builtin_popcountll(tree[t_pos - 1]);
  if(t_pos == b_tree) {
    if(c_root > 1)
      return z->v[z_pos] & TSIZEMASK;
    else
      return t_size - 1;
  }

  // search for wich tree you are
  uint64_t z_m = c_root - 1;
  uint64_t t_m = 0;
  for(uint64_t i = z_pos; i < z_pos + c_root - 1; i++) {
    uint64_t st_size = z->v[i] & TSIZEMASK;
    uint64_t z_skip = z->v[i] >> BITSxTSIZE;

    // found the tree
    if(t_pos + t_m + st_size > b_tree) {
      // move to the first child of this subtree
      return get_size_rec(tree, z, t_pos + t_m + 1, st_size, z_pos + z_m, b_tree);
    }

    // is the brother but not the last
    if(t_pos + t_m + st_size == b_tree && i < z_pos + c_root - 2) {
      return z->v[i + 1] & TSIZEMASK;
    // is the last
    } else if(t_pos + t_m + st_size == b_tree) {
      return t_size - t_m - st_size - 1;
    }
    // moving to next tree
    t_m += st_size;
    z_m += z_skip;
  }
  // is in the last child of a node
  return get_size_rec(tree, z, t_pos + t_m + 1, t_size - t_m - 1, z_pos + z_m, b_tree);
}

static uint64_t get_size(uint8_t *tree, uint64_t t_size, vu64_t *z, uint64_t b_tree) {
  if(b_tree == 0) return t_size;

  uint64_t t_pos = 1;
  uint64_t z_pos = 0;

  return get_size_rec(tree, z, t_pos, t_size, z_pos, b_tree);
}

node_t k2read_node__(const k2mat_t *m, size_t p)
{
  p+=m->offset;
  assert(p<m->pos && p<m->lenb);
  if(p%2==0)
    return m->b[p/2] & 0xF;
  else 
    return  (m->b[p/2] >> 4) & 0xF;
}

size_t k2add_node__(k2mat_t *m, node_t n)
{
  assert(!m->read_only);
  assert(n<ILLEGAL_NODE);
  assert(m->lenb%2==0);            // #positions must be even
  // make sure there is space
  if(m->pos >= m->lenb) {
    assert(m->pos ==m->lenb);
    m->lenb = 16+2*m->lenb;          // more than double number of positions 
    m->b = (uint8_t*) realloc(m->b, m->lenb/2); // each byte stores two positions
    //if(m->b==NULL) quit("Unable to enlarge k2-tree",__LINE__,__FILE__);
  }
  assert(m->pos<m->lenb);
  // since a node is stored in 4 bits, we store two nodes in a byte
  if(m->pos%2==0)
    m->b[m->pos/2] = (uint8_t) n;    // note: we are writing 0 at m->pos+1, it's ok we are at the end 
  else
    m->b[m->pos/2] = (m->b[m->pos/2] & 0x0F) | (uint8_t) (n<<4);
  // return position where node was stored and advance by 1
  return m->pos++; 
}

k2mat_t compress_k2mat_t(size_t size, size_t asize, k2mat_t* a,
                         uint32_t* P_size, uint32_t** P, // store pointers
                         uint32_t* rank_size, uint32_t** rank_p, uint32_t rank_block) { // store rank 0000

  uint64_t lvs = ceil_log2((uint64_t)asize);

  vu64_t z;
  vu64_init(&z);
  size_t pos = 0;

  k2dfs_sizes(asize, a, &pos, &z, (uint32_t) lvs);

  uint8_t *text = (uint8_t*) malloc(sizeof(uint8_t) * a->pos);

  for(size_t i = 0; i < a->pos; i++) {
    text[i] = (uint8_t) k2read_node__(a, i);
  }

  int64_t *csa = malloc(sizeof(uint64_t) * a->pos);
  int64_t *plcp = malloc(sizeof(uint64_t) * a->pos);
  int64_t *lcp = malloc(sizeof(uint64_t) * a->pos);

  if(libsais64(text, csa, a->pos, 0, NULL) != 0) {}
  if(libsais64_plcp(text, csa, plcp, a->pos) != 0) {}
  if(libsais64_lcp(plcp, csa, lcp, a->pos) != 0) {}

  dsu u;
  dsu_init(&u, a->pos);

  st_t st;
  st_init(&st, a->pos, lcp);

  size_t* prev = malloc(sizeof(size_t) * (a->pos + 1));
  for(size_t i = 0; i < (a->pos + 1); i++) prev[i] = -1;

  // information variables
  uint64_t amount_of_groups = 0;

  for(size_t i = 0; i < a->pos; i++) {
    uint64_t curr_start_pos = csa[i];
    uint64_t curr_end_pos = csa[i] + get_size(text, a->pos, &z, csa[i]) - 1;
    if(prev[curr_end_pos - curr_start_pos + 1] == -1) {
      prev[curr_end_pos - curr_start_pos + 1] = i;
      continue;
    }

    int64_t pos = prev[curr_end_pos - curr_start_pos + 1];
    uint64_t prev_start_pos = csa[pos];

    // ignoring leaves
    if(curr_end_pos - curr_start_pos + 1 <= 32 / 4) {
      prev[curr_end_pos - curr_start_pos + 1] = i;
      continue;
    }

    // check that the tree are same length
    if(st_query(&st, pos + 1, i) < curr_end_pos - curr_start_pos + 1) {
      amount_of_groups++;
    } else {
      dsu_union_set(&u, curr_start_pos, prev_start_pos);
    }
    prev[curr_end_pos - curr_start_pos + 1] = i;
  }

  free(csa);
  free(plcp);
  free(lcp);
  free(prev);
  st_free(&st);

  k2mat_t ca = K2MAT_INITIALIZER;

  uint64_t* prefix_help = (uint64_t*) malloc(sizeof(uint64_t) * a->pos);
  for(size_t i = 0; i < a->pos; i++) prefix_help[i] = 0;


  vu64_t P_h;
  vu64_init(&P_h);

  for(size_t i = 0; i < a->pos; i++) {
    uint64_t repre = dsu_find_set(&u, i);
    if(repre != i) { // cancel identical subtree
      vu64_grow(&P_h, 1);
      P_h.v[P_h.n - 1] = (uint32_t) repre - prefix_help[repre - 1];
      k2add_node__(&ca, 0);
      size_t next_i = i + get_size(text, a->pos, &z, i) - 1;
      for(size_t fill = i; fill <= next_i; fill++) {
        prefix_help[fill] = prefix_help[fill - 1] + 1;
      }
      prefix_help[next_i]--;
      i = next_i;
    } else {
      if(i > 0) prefix_help[i] = prefix_help[i - 1];
      k2add_node__(&ca, text[i]);
    }
  }
  
  free(prefix_help);

  *P = (uint32_t*) malloc(sizeof(uint32_t) * P_h.n);
  *P_size = (uint32_t) P_h.n;
  for(size_t i = 0; i < P_h.n; i++) {
    assert(ca.pos > P_h.v[i]);
    (*P)[i] = (uint32_t) P_h.v[i];
  }

  vu64_free(&P_h);
  vu64_free(&z);
  free(text);
  dsu_free(&u);

  *rank_size = (uint32_t) (ca.pos + rank_block - 1) / rank_block + 1;
  *rank_p = (uint32_t*) malloc(sizeof(uint32_t) * (*rank_size));
  uint32_t sum = 0;
  for(size_t i = 0; i < ca.pos; i++) {
    if(i % rank_block == 0) {// finish block
      (*rank_p)[(i + rank_block - 1) / rank_block] = sum;
    }
    sum += k2read_node__(&ca, i) == 0;
  }

  (*rank_p)[(ca.pos + rank_block - 1) / rank_block] = sum;

  return ca;
}


typedef uint64_t minimat_t; // explict matrix aka minimat (currently 2x2)

static uint32_t MMsize = 0;  // a certainly incorrect value 
static uint32_t Minimat_node_ratio = 0;  // certainly incorrect value

// multiplication functions for minimats of different sizes

// global variable containing the product of every pair of possible minimatrices
// to be initialized in minimats_init();
static minimat_t k2read_minimat__(const k2mat_t *b, size_t *p) {
  if(MMsize==2) {
    assert(Minimat_node_ratio==1);
    return (minimat_t) k2read_node__(b,(*p)++);
  }
  else {
    assert(MMsize==4);
    assert(Minimat_node_ratio==4);
    minimat_t res = 0;
    for(int i=0;i<Minimat_node_ratio;i++) {
      res <<= 4;
      res |= (minimat_t) k2read_node__(b,(*p)++);
    }
    return res;
  }    
}

uint32_t rank_p(const k2mat_t *m, size_t pos, uint32_t R_size, uint32_t block_size, uint32_t *rank_h) {
  uint32_t block = (uint32_t) pos / block_size;

  assert(block < R_size);
  uint32_t ret = rank_h[block];

  for(uint32_t to_read = block * block_size; to_read < pos; to_read++) {
    ret += k2read_node__(m, to_read) == 0;
  }
  return ret;
}

void k2dfs_visit__(size_t size, const k2mat_t *m, size_t *pos, size_t *nodes, size_t *minimats, size_t *nz, size_t *onodes, size_t *ominimats,
                   uint32_t P_size, uint32_t *P, uint32_t R_size, uint32_t block_size, uint32_t *rank_h)
{
  assert(size>MMsize);
  assert(size%2==0);
  assert(*pos<m->pos); // implies m is non-empty
  node_t root = k2read_node__(m,*pos); (*pos)++;
  (*nodes)++;
  (*onodes)++;
  if(root==ALL_ONES) { // just coincidence that ALL_ONES == POINTER
    (*onodes)--;
    uint32_t aux = (uint32_t) *pos; // remember where to comback
    size_t aux_nodes = *nodes; // to not overcount nodes
    size_t aux_minimats = *minimats; // to not overcount minimats

    uint32_t rp = rank_p(m, *pos - 1, R_size, block_size, rank_h);
    assert(rp < P_size);
    *pos = P[rp];
    assert(*pos < m->pos);
    k2dfs_visit__(size, m, pos, nodes, minimats, nz, onodes, ominimats, P_size, P, R_size, block_size, rank_h); // read submatrix and advance pos
    
    *nodes = aux_nodes; // get back correct stats
    *minimats = aux_minimats; // get back correct stats
    
    // moving back to pointer
    *pos = aux;
    return; // all 1's matrix consists of root only
  }
  for(int i=0;i<4;i++) 
    if(root & (1<<i)) {
      if(size==2*MMsize) { // end of recursion
        (*minimats)++;
        (*ominimats)++;
        minimat_t mm = k2read_minimat__(m,pos); // read minimat and advance pos
        *nz += (size_t) __builtin_popcountll(mm);
      }
      else { // recurse on submatrix
        k2dfs_visit__(size/2, m, pos, nodes, minimats, nz, onodes, ominimats, P_size, P, R_size, block_size, rank_h); // read submatrix and advance pos
      }
    }
}

bool k2is_empty__(const k2mat_t *m)
{
  return m->pos == 0;
}

static int mstats__(size_t asize, const k2mat_t *a, size_t *pos, size_t *nodes, size_t *minimats, size_t *nz, size_t *onodes, size_t *ominimats,
                    uint32_t P_size, uint32_t *P, uint32_t R_size, uint32_t block_size, uint32_t *rank_h)
{
  *pos=*nodes=*minimats=*nz=*onodes=*ominimats=0;
  if(!k2is_empty__(a)) k2dfs_visit__(asize,a,pos,nodes,minimats,nz, onodes, ominimats, P_size, P, R_size, block_size, rank_h);
  int eq = mequals(asize,a,a);
  assert(eq<0);
  return -eq;
}

// write to :file statistics for a k2 matrix :a with an arbitrary :name as identifier
size_t mshow_stats__(size_t size, size_t asize, const k2mat_t *a, const char *mname,FILE *file,
                     uint32_t P_size, uint32_t *P, uint32_t R_size, uint32_t block_size, uint32_t *rank_h) {
  size_t pos, nodes, minimats, nz, onodes, ominimats;
  fprintf(file,"%s:\n matrix size: %zu, leaf size: %d, k2_internal_size: %zu\n",mname,size,MMsize,asize);  
  int levels = mstats__(asize,a,&pos,&nodes,&minimats,&nz, &onodes, &ominimats, P_size, P, R_size, block_size, rank_h);
  assert(pos==nodes+minimats*Minimat_node_ratio); // check that the number of positions is correct
  fprintf(file," Levels: %d, Nodes: %zu, Minimats: %zu, Nonzeros: %zu, Original Nodes: %zu, Original Minimats: %zu\n",
          levels,nodes,minimats, nz, onodes, ominimats);
  // each pos takes 4 bits, so three size in bytes is (pos+1)/2         
  fprintf(file," Tree size: %zu bytes, Bits x nonzero: %lf\n",
          (pos+1)/2 , 4.0*(double)(pos)/(double) nz);
  return nz;
}

void print_ck2(k2mat_t *a, uint32_t P_size, uint32_t *P, uint32_t R_size, uint32_t *rank_h) {
  for(size_t i = 0; i < a->pos; i++) {
    for(size_t bit = 0; bit < 4; bit++) {
      if(i % 2)
        printf("%d", (((a->b[i / 2] >> 4) & (1 << bit)) > 0));
      else
        printf("%d", (((a->b[i / 2] & 15) & (1 << bit)) > 0));
    }
    printf(" ");
  }
  printf("\n");
  for(size_t i = 0; i < P_size; i++) {
    printf("%" PRIu32 " ", P[i]);
  }
  printf("\n");
  for(size_t i = 0; i < R_size; i++) {
    printf("%" PRIu32 " ", rank_h[i]);
  }
  printf("\n");
}

void check_rank_info(k2mat_t *ca, uint32_t R_size, uint32_t block_size, uint32_t *rank_h) {
  int prefix_sum = 0;
  for(size_t i = 0; i < ca->pos; i++) {
    uint8_t node = (uint8_t) k2read_node__(ca, i);
    assert(rank_p(ca, i, R_size, block_size, rank_h) == prefix_sum);
    prefix_sum += node == 0;
  }
}

int main(int argc, char* argv[]) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  k2mat_t a = K2MAT_INITIALIZER;
  size_t size, asize, totnz=0;
  uint32_t rank_block = 64;
  uint32_t P_size;
  uint32_t* P = NULL;
  uint32_t rank_size;
  uint32_t* rank_p = NULL; 

  MMsize = 2;
  Minimat_node_ratio = 1;
  char* k2name_file = NULL;
  int c;
  while((c = getopt(argc, argv, "b:f:")) != -1) {
    switch (c) {
      case 'b':
        rank_block = atoi(optarg); break;
      case 'f':
        k2name_file = optarg; break;
      case '?':
        fprintf(stderr, "Unknown option %c\n", optopt);
        exit(1);
    }
  }

  if(k2name_file == NULL) {
    fprintf(stderr, "Missing file argument: -f <name file>\n");
    exit(1);
  }

  size = mload_from_file(&asize, &a, k2name_file); // also init k2 library
  totnz += mshow_stats(size, asize, &a, basename(k2name_file), stdout);


  k2mat_t ca = compress_k2mat_t(size, asize, &a, &P_size, &P, &rank_size, &rank_p, rank_block);
  check_rank_info(&ca, rank_size, rank_block, rank_p);
  assert(P_size == rank_p[(ca.pos + rank_block - 1) / rank_block]);

  uint64_t bits_ca = sizeof(ca) + sizeof(ca.lenb) + sizeof(ca.offset)
                 + sizeof(ca.pos) + sizeof(ca.read_only) + sizeof(int8_t) * ((ca.pos + 1) / 2);
  bits_ca *= 8;
  uint64_t bits_rank = sizeof(uint32_t) * rank_size + 1;
  bits_rank *= 8;
  uint64_t bits_P = sizeof(uint64_t) * P_size + 1;
  bits_P *= 8;
  uint64_t bits_extra = bits_rank + bits_P;

  fprintf(stdout, "COMPRESSED TREE\n");
  fprintf(stdout, "Nodes: %ld\n", ca.pos);
  fprintf(stdout, "COMPRESSED SPACE\n");
  fprintf(stdout, "Bits compressed k2pdf: %" PRIu64 " Removed Trees: %" PRIu32 "\nBits extra info: %" PRIu64 " bits rank: %" PRIu64 " bits P: %" PRIu64
                  "\nTotal bits: %" PRIu64 " Total bytes: %" PRIu64 "\n",
                  bits_ca, P_size, bits_extra, bits_rank, bits_P, bits_ca + bits_extra, (bits_ca + bits_extra) / 8);

  //print_ck2(&a, 0, NULL, 0, NULL);
  //print_ck2(&ca, P_size, P, rank_size, rank_p);
  mshow_stats__(size, asize, &ca, basename(k2name_file), stdout, P_size, P, rank_size, rank_block, rank_p);

  char file_ck2[strlen(k2name_file) + 5];
  strcpy(file_ck2, k2name_file);
  file_ck2[strlen(k2name_file) - 2] = 'c';
  file_ck2[strlen(k2name_file) - 1] = 'k';
  file_ck2[strlen(k2name_file)] = '2';
  file_ck2[strlen(k2name_file) + 1] = '\0';
  msave_to_file(size, asize, &ca, file_ck2);

  char file_p[strlen(file_ck2) + 4];
  strcpy(file_p, file_ck2);
  file_p[strlen(file_ck2)] = '.';
  file_p[strlen(file_ck2) + 1] = 'p';
  file_p[strlen(file_ck2) + 2] = '\0';
  FILE *fp = fopen(file_p, "w");  
  fwrite(&P_size, sizeof(uint32_t), 1, fp);
  fwrite(P, sizeof(uint64_t), P_size, fp);
  fclose(fp);

  char file_r[strlen(file_ck2) + 4];
  strcpy(file_r, file_ck2);
  file_r[strlen(file_ck2)] = '.';
  file_r[strlen(file_ck2) + 1] = 'r';
  file_r[strlen(file_ck2) + 2] = '\0';
  FILE *fr = fopen(file_r, "w");  
  fwrite(&rank_size, sizeof(uint32_t), 1, fr);
  fwrite(&rank_block, sizeof(uint32_t), 1, fr);
  fwrite(rank_p, sizeof(uint32_t), P_size, fr);
  fclose(fr);

  free(P);
  free(rank_p);
  matrix_free(&ca);
  matrix_free(&a);
  minimat_reset();
  return 0;
}
