# k2tree binary matrix representation


This repository contains a set of functions for working with square sparse binary matrices represented with a k$^2$-tree. Matrix sum and product operations use boolean algebra with the scalar operations `+` and `*` replaced by logical `or` and `and` respectively. A version supporting $Z_2$ arithmetic can be easily derived if needed. 


## Prerequisites 

A modern `gcc` supporting `c11` and `make`



## Installation 

Clone/download the repostory then `make release`



## Getting started


### The uncompressed matrix format


The input matrix must be square and represented using one byte per entry (so the uncompressed file lenght is the square of the matrix dimension). Each entry should be a 0 or 1 value (ie `\x00` or `\x01`). The file `t8.bbm` contains an 8x8 matrix in this format, type `od -td1 -w8 -An -v t8.bbm` to see its content:
```
    0    0    0    0    0    1    1    1
    0    0    0    0    0    0    1    0
    0    0    0    0    0    0    1    1
    0    0    0    0    0    1    0    1
    1    1    1    1    0    0    1    0
    1    1    1    1    0    0    0    1
    1    1    1    1    0    0    0    1
    1    1    1    1    0    0    1    0
```
The extension `.bbm` stands for **b**inary **b**yte **m**atrix.



### Conversion to/from k2tree representation 

The executable `k2comp.x` is used to compress and decode a single matrix to/from a k2tree representation; run it without arguments to get basic usage instructions
```
Usage:
	  k2comp.x [options] filename

Options:
	-d      decompress
	-n      do not write output file, only show stats
	-m M    minimatrix size (def. 2), compression only
	-e ext  outfile extension (def. compr: ".k2", decompr: ".d")
	-c      compress->decompress->check
	-v      verbose

Default action is to compress filename to filename.k2
``` 


Invoked without the `-d` or `-c` switch `k2comp.x` compresses the input matrix to the k2 format. 
The input file must be a binary square matrix represented with one byte per entry (so the file size must be a perfect square). 

When invoked with `-d` the input file must be a k2 matrix which is then expanded into the one byte per entry format. 

When invoked with `-c` the program compresses the input matrix, then decompresses it and verify that the decompressed matrix matches the original matrix. 


The option `-m` can be used only in compression, currently only with one of the two values `2` and `4`. This parameter is the size of the matrices stored at the leaf of the k2 tree (except the leaves representing the submatrices of all 1's which can be of any size). A large leaf size usually yields larger files but improves the running time for the aritmetic operations over the matrices (as the tree is shallower). The value of the parameter `-m` is stored in the k2 format so it does not have to be provided for decompression.


### Product of matrices in k2 format

The executable `k2mult.c` can be used to multiply two compressed matrices in k2 format. The matrices must have the same size and must have been compressed with the same `-m` parameter. The output is still in k2 format.

```
Usage:
	  k2mult.x [options] iname1 iname2

Options:
	-n      do not write output file, only show stats
	-e ext  extension for the output file (def. .prod)
	-c      check multiplication (can be slow for large matrices!)
	-v      verbose

```
When invoked with `-c`, after computing the product in k2 format, the program uncompresses the input matrices and the product and verify that the uncompressed product is identical to the product computed with the tradizional algorithm applied to the uncompressed inputs. For large matrice this verification can be slow and space consuming.   


### Example:

The following sequence of compression, matrix multiplication, decompression and display operations
```
k2comp.x t8.bbm
k2mult.x t8.bbm.k2 t8.bbm.k2
k2comp.x -d t8.bbm.k2.prod
od -td1 -w8 -An -v t8.bbm.k2.prod.d
```
should eventually display the matrix `t8.bbm` squared:
```
    1    1    1    1    0    0    1    1
    1    1    1    1    0    0    0    1
    1    1    1    1    0    0    1    1
    1    1    1    1    0    0    1    1
    1    1    1    1    0    1    1    1
    1    1    1    1    0    1    1    1
    1    1    1    1    0    1    1    1
    1    1    1    1    0    1    1    1
```


### Matrices represented as bitarrays

The library contains the code also for compressing and operating binary matrices using a bitarray, ie using one bit per entry plus a small overhead. To make the interchange between the two compressed formats very simple, the callable functions (whose prototypes are in `k2.h` and `b128.h`) have the same names. Hence, a program using the k2 format can be transformed into one using the bitarray format by redefining a few constats. See the use of the `B128MAT` compilation constant in the source files `k2comp.c` and `k2mult.c` and in the `makefile`.

The programs `b128comp.x` and `b128mult.x` work exactly like  `k2comp.c` and `k2mult.c` except that they use the bitarray representation instead of the k2 format. 


### Product of uncompressed matrices

The program `bbmmult.x` computes the product of two uncompressed matrices using `openmp` to speedup the computation. This tool has been provided mainly as faster alternative to the `-c` option to check the correctness of `k2mult.x` and `b128mult.x`.



### Additional tools 

The files `k2stat.sh`, `k2test.sh`, `k2square.sh` are bash script designed to test `k2comp.x` and `k2mult.x` on a set of input files.  
