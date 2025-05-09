# k2tree binary matrix representation


This repository contains a set of functions for working with square sparse boolean matrices represented with a $k^2$-tree. Matrix sum and product operations use boolean algebra with the scalar operations `+` and `*` replaced by logical `or` and `and` respectively. A version supporting $Z_2$ arithmetic can be easily derived if needed. 


## Prerequisites 

A modern `gcc` supporting `c11` and `make`



## Installation 

Clone/download the repostory then `make release`

All tools invoked without arguments provide basic usage instructions. 



## Getting started


### The uncompressed matrix formats

The simplest uncompressed format is the **b**inary **b**yte **m**atrix (extension `.bbm`). In this format the input matrix must be square and represented using one byte per entry (so the `.bbm` file lenght is the square of the matrix dimension). Each entry should be a 0 or 1 value (ie `\x00` or `\x01`). The file `t8.bbm` contains an 8x8 matrix in this format, type `od -td1 -w8 -An -v t8.bbm` to see its content:
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

The advantage of this format is that we can look inside a file with command line tools (like `od`), but it is extremely space inefficient since it does not take advantage of sparsity an uses one byte for each 0/1 entry.

Another popular uncompressed format consists of a text file containing only the positions of the nonzero entries. Each line should contain the row and column indexes of a single entry written in decimal and separated by a whitespace character.  The same entry should not appear twice, but the order of entries can be arbitrary. Indexes are 0-based so the matrix above could be represented by the text file `t8.txt`:
```
0 5
0 6
0 7
1 6
3 5
2 6
2 7
3 7
4 0
4 1
4 2
4 3
5 0
5 1
5 2
5 3
6 0
6 1
6 2
6 3
7 0
7 1
7 2
7 3
4 6
5 7
6 7
7 6
```

### Conversion to/from .bbm representation 

The executable `k2bbm.x` is used to compress and decode a single `.bbm` matrix to/from a k2tree representation; run it without arguments to get basic usage instructions
```
Usage:
      k2bbm.x [options] filename

Options:
    -d      decompress
    -n      do not write the output file, only show stats
    -o out  outfile name (def. compr: infile.k2, decompr: outfile.bbm)
    -m M    minimatrix size (def. 2), compression only
    -1      do not compact all 1's submatrices, compression only
    -c      compress->decompress->check
    -h      show this help message
    -v      verbose

Default action is to compress filename to filename.k2
``` 


Invoked without the `-d` or `-c` switch `k2bbm.x` compresses the input `.bbm` matrix to the k2 format. 

When invoked with `-d` the input file must be a k2 matrix which is then expanded into the `.bbm` format. 

When invoked with `-c` the program compresses the input matrix, then decompresses it and verify that the decompressed matrix matches the original matrix. 

The option `-m` can be used only in compression, currently only with one of the two values `2` and `4`. This parameter is the size of the matrices stored at the leaf of the k2 tree (except the leaves representing the submatrices of all 1's which can be of any size). A large leaf size usually yields larger files but improves the running time for the aritmetic operations over the matrices (as the tree is shallower). The value of the parameter `-m` is stored in the k2 format so it does not have to be provided for decompression.


### Conversion to/from textual one line per entry sparse representation 

The executable `k2sparse.x` is used to compress and decode a single textual matrix to/from a k2tree representation; run it without arguments to get basic usage instructions
```
Usage:
      k2sparse.x [options] filename

Options:
	-d      decompress
	-n      do not write the output file, only show stats
	-o out  outfile name (def. compr: infile.k2, decompr: infile.txt)
	-s S    matrix actual size (def. largest index+1) [compression only]
	-m M    minimatrix size (def. 2) [compression only]
	-1      compact all 1's submatrices [compression only]
	-c      compress->decompress->check
	-h      show this help message
	-v      verbose

Default action is to compress filename to filename.k2
```

Invoked without the `-d` or `-c` switch `k2sparse.x` compresses a textual matrix to the k2 format. 

When invoked with `-d` the input file must be a k2 matrix which is then expanded into the textual format. 

When invoked with `-c` the program compresses the input matrix, then decompresses it and verify that the decompressed matrix matches the original matrix. Since the order of the entries in the textual file is arbitrary the verification involves sorting and searching and is done invoking the tool `matrixcmp.x`.

The option `-m` can be used only in compression, currently only with one of the two values `2` or `4`. This parameter is the size of the matrices stored at the leaf of the k2 tree (except the leaves representing the submatrices of all 1's which can be of any size). A large leaf size usually yields larger files but improves the running time for the aritmetic operations over the matrices (as the tree is shallower). The value of the parameter `-m` is stored in the k2 format so it does not have to be provided for decompression.

By default the input matrix is assumed to be of size 1+(largest index in the input file). The option `-s` can be used to force the size of the input matrix to a specific (larger) value. Because of the algorithm used to compress textual matrices, currently the largest admissible matrix size is $2^{32}$; this limitation can be removed if needed using a slightly more complex compression algorithm. 



## Product of matrices in k2 format

The executable `k2mult.c` can be used to multiply two compressed matrices in k2 format. The matrices must have the same size and must have been compressed with the same `-m` parameter. The output is still in k2 format.

```
Usage:
      k2mult.x [options] infile1 infile2

Options:
    -n        do not write output file, only show stats
    -o out    outfile name (def. infile1.prod)
    -1        compact all 1's submatrices in the result matrix
    -i info   infile1 subtree info file
    -j info   infile2 subtree info file
    -c        check multiplication (O(n^3) time and O(n^2) space!)
    -h        show this help message
    -v        verbose

Multiply two compressed matrices stored in infile1 and infile2
```
When invoked with `-c`, after computing the product in k2 format, the program uncompresses the input matrices and the product and verify that the uncompressed product is identical to the product computed with the traditional $O(n^3)$ time algorithm applied to the uncompressed inputs. For large matrices this verification can be slow and space consuming. 
For an explanation of the `-i` and `-j` options see section *Enriched format* below. 


### Example:

The following sequence of compression, matrix multiplication, decompression and display operations
```
k2bbm.x t8.bbm
k2mult.x t8.bbm.k2 t8.bbm.k2
k2bbm.x -d t8.bbm.k2.prod
od -td1 -w8 -An -v t8.bbm.k2.prod.bbm
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


## Enriched format

In order to speedup operations on compressed matrices in k2 format it is possible to use some extra information 
on its largest subtrees. This information must be computed with `k2subtinfo.x`, for example:
```bash
k2subtinfo.x -vc m.k2 -o m.k2.info 
```
computes the information for the compressed matrix `m.k2` and stores it in `m.k2.info`. By default the information is computed for the subtrees whose size is at least the square root of the size of the k2 tree representing the whole matrix. 
Use the option `-N limit` to compute the info only for subtrees of size larger than `limit`, or use the option `-D d` to compute the info only for the subtree at depth up to `d`. 

At the moment this additional information can be used only to speed-up matrix-matrix multiplication. For example write
```bash
k2mult.x m1.k2 m2.k2 -i m1.k2.info -j m2.k2.info 
```
to compute the product `m1.k2` times `m2.k2` using the extra information `m1.k2.info` for the compressed matrix `m1.k2`
and `m2.k2.info` for the compressed matrix `m2.k2` (the information can be used also for only one of the two input matrices).

NOTE: to see a real speed improvement at the moment it is necessary to use the release version (`make release`). 



## Pagerank computation 

The Pagerank computation is ideal for testing the speed of the matrix-vector product in a real-world scenario.
Assuming that the input (web) matrix is given in mtx format it is first necessary to preprocess 
it using the `mtx2rowm` tool that after the conversion provides some minimal instructions to
compress the input matrix (using `k2sparse.x` or `k2blockc.py`) and later compute the Pagerank
vector (using `k2pagerank.x`) possibly using multiple threads. (The program `k2bpagerank.x` is an experiment
where thread syncronization is done using pthread barriers, but the performances are very similar). 



## Matrices represented as bitarrays

The library also contains the code for compressing and operating on boolean matrices using a bitarray, ie using one bit per entry plus a small overhead. To make the conversion between the two compressed formats very simple, the callable functions (whose prototypes are in `k2.h` and `b128.h`) have the same names. Hence, a program using the k2 format can be transformed into one using the bitarray format by redefining a few constants. See the use of the `B128MAT` compilation constant in the source files `k2bbm.c` and `k2mult.c` and in the `makefile`. Creation of bitarray matrices is currently not supported for textual input matrices. Since bitarray representation does not take advantage of sparsity, the largest supported size is $2^{30}$. 

The programs `b128sparse.x`, `b128bbm.x` and `b128mult.x` work exactly like  `k2sparse.x`, `k2bbm.x` and `k2mult.x` except that they use the bitarray representation instead of the k2 format. 


## Product of bbm matrices with openmp

The program `bbmmult.x` computes the product of two `.bbm` matrices using `openmp` to speedup the computation. This tool has been provided mainly as an alternative to the `-c` option to check the correctness of `k2mult.x` and `b128mult.x`.



## Additional tools 


The program `k2info.x` display statics on the k2-compressed files passed on the command line.

The script `submatrix.py` can be used to extract a square submatrix form a matrix in textual form.

The files `k2test.sh`, `k2btest.sh`, `k2square.sh` are bash script designed to test `k2sparse.x`, `k2bbm.x` and `k2mult.x` on a set of input files.  

