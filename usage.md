# Compression formats and operations 


## Plain Depth Dirst 

Plain Depth First representation with compression of all 1 submatrices

### compression

```
./k2sparse.x m.txt -o m.k2
```

### squaring

```
./k2mult.x -q m.k2 -o msquare.k2
```

### multiplication
```
./k2mult.x m1.k2 m2.k2 -o m3.k2
```


## Enriched Depth First 

As above but uses subtree size information for the upper levels

### compression

```
./k2sparse.x m.txt -o m.k2
./k2subtinfo.x m.k2 
```
see `./k2subtinfo.x -h` for the additional options `-D`, `-N`, `-M` to control the generation of subtree information. 


### squaring

```
./k2mult.x -q m.k2 -i m.k2.sinfo -o msquare.k2
```

### multiplication
```
./k2mult.x m1.k2 m2.k2 -i m1.k2.sinfo -j m2.k2.sinfo -o m3.k2
```

### squaring/multiplication with dynamic subtree info

To compute the subtree info on the fly for small subtrees use option `-e`; for example
```
./k2mult.x -e -q m.k2 -i m.k2.sinfo -o msquare.k2
```


## Compressed Depth First 

Compress the $k^2$ tree by detecting identical subtrees. No automatic compression of submatrices of all 1s.

###  Compression

```
./k2sparse.x m.txt -o mx.k2
./k2cpdf.x m.k2
./k2subtinfo.x m.ck2 
```


### squaring

```
./k2mult.x -q m.ck2 -i m.ck2.sinfo -I m.ck2.p -o msquare.k2
```

### multiplication
```
./k2mult.x m1.k2 m2.k2 -i m1.k2.sinfo -I m1.ck2.p -j m2.k2.sinfo -J m2.ck2.p -o m3.k2
```

### squaring/multiplication with dynamic subtree info

As above, to compute the subtree info on the fly for small subtrees use option `-e`; for example
```
./k2mult.x -e -q m.ck2 -i m.ck2.sinfo -I m.ck2.p -o msquare.k2
```


# Available C functions 

These are the main functions handling k2matrices, see the above tools for example of usage

```C
// check if two matrices are equal element by element
bool mequals(const k2mat_t *a, const k2mat_t *b);
// creates a size x size zero matrix
k2mat_t mat_zero(const k2mat_t *a);
//creates a size x size identity matrix
k2mat_t mat_identity(const k2mat_t *a);
// add indentity matrix to a
void madd_identity(k2mat_t *a);
// sum (logical or) of two k2 matrices a and b writing the result to c
void msum(const k2mat_t *a, const k2mat_t *b, k2mat_t *c);
// multiply two k2 matrices a and b writing the result to c
void mmult(const k2mat_t *a, const k2mat_t *b, k2mat_t *c);
// right multiply a k2 matrix :a by a vector :x writing the result to :y
void mvmult(const k2mat_t *a, double *x, double *y, bool clear_y);
// read/write a matrix in text format
size_t mread_from_textfile(k2mat_t *a, char *iname, size_t xsize);
void mwrite_to_textfile(const k2mat_t *a, char *outname);
// load/save a k2-matrix in compressed form
size_t mload_from_file(k2mat_t *a, const char *filename);
void msave_to_file(const k2mat_t *a, const char *filename);
```

