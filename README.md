# k2tree binary matrix representation


This repository contains a library for working with square sparse binary matrices. The matrix elements are considered logical values so the operations `+` and `*` are repalced by the logival `or` and `and`  respectively. 


## Prerequisites 

A modern `gcc` and `make`



## Installation 

Clone/download the repostory then `make release`



## Getting started


### The uncompressed matrix format


The input matrix must be square and represented using one byte per entry (so the uncompressed file lenght is the square of the matrix dimension). Each entry should be a 0 or 1 value (ie `\x00` or `\x01`). The file `t8.bin` contains an 8x8 matrix in this format, type `od -td1 -w8 -An -v t8.bin` to see its content:
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



### Converion to/from k2tree representation 


The executable `k2comp.x` is used to compression and decode a single matrix to/from a k2tree representation, run them without arguments to get basic usage instructions
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

When invoked with `-d` the input file must be a k2 matrix which is then espanded into the one byte per entry format. 

When invoked with `-c` the program compresses the input matrix, then decompresses it and verify that the decompressed matrix matches the original matrix. 


The option `-m` can be used only in compression, currently only with one of the two values `2` and `4`. This parameter is the size of the matrices stored at the leaf of the k2 tree (except the leaves representing the submatrices of all 1's which can be of any size). A large leaf size usually yields larger files but improves the running time for the aritmetic operations over the matrices (as the tree is shallower). The valu of the parameter `-m` is stored in the k2 format so it does not have to be provided for decompression.


### Product of matrices in k2 format

The executable `k2mult.c` is used to multiply two compressed matrices in k2 format. The matrices must have the same size and must have been compressed with the same `-m` parameter. The output is still in k2 format.

```
Usage:
	  k2mult.x [options] iname1 iname2

Options:
	-n      do not write output file, only show stats
	-e ext  extension for the output file (def. .prod)
	-c      check multiplication (can be slow for large matrices!)
	-v      verbose

```
When invoked with `-c` after computing the product in k2 format, the program uncompressed the input matrices and the product and verify that the uncompressed products is identical to the product computed with the tradizional algorithm applied to the uncompressed inputs. For large matrice this verification can be slow and space consuming.   


### Example:

The following sequence of compression, matrix multiplication, decompression and display
```
k2comp.x t8.bin
k2mult.x t8.bin.k2 t8.bin.k2
k2comp.x -d t8.bin.k2.prod
od -td1 -w8 -An -v t8.bin.k2.prod.d
```
should display the matrix `t8.bin` squared
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

The library contains the code also for compressing and operating binary matrices using a bitarray, ie using one bit per entry plus a small overhead. To make the interchange between the two compressed formats very simple the external functions (whose prototypes are in `k2.h` and `b128.h`) have the same names. Hence, a program using the k2 format can be transormed into one using the bitarray format by redefining a few constats. See the use of the `B128MAT` compilation constant in the source files `k2comp.c` and `k2mult.c` and in the `makefile`.

The programs `b128comp.x` and `b128mult.x` works exactly like  `k2comp.c` and `k2mult.c` except that the use the bitarray representation instead that the k2 format. 


### Product of uncompressed matrices

The program `bbmmult.x` computes the product of two uncompressed matrices using `openmp` to speedup the computation. This tool has been provided mainly as faster alternative to the `-c` option to check the correctness of `k2mult.x` and `b128mult.x`.



### Additional tools 

The files `k2stat.sh`, `k2test.sh`, `k2square.sh` are bash script designed to test `k2comp.x` and `k2mult.x` on a set of input files.  
