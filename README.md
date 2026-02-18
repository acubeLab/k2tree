# k2tree binary matrix representation


This repository contains a set of functions for working with square sparse boolean matrices represented with a $k^2$-tree. Matrix sum and product operations use boolean algebra with the scalar operations `+` and `*` replaced by logical `or` and `and` respectively. A version supporting $Z_2$ arithmetic can be easily derived if needed. 


## Prerequisites 

A modern `gcc` supporting `c11` and `make`

A cmake version `>= 3.10`

## Installation 

Clone/download the repostory then:

```
cd libsais; cmake .; make; cd ..
make release
```

All tools invoked without arguments provide basic usage instructions. 

## Getting started

### TL;DR

Given the text file `m.txt` containing in each line the row and column index of a nonzero element (0-based):
```
k2sparse.x m.txt -o m.k2
```
compute $k^2$-tree representation of `m.txt`.
```
k2mult.x m.k2 m.k2 -o m2.k2
```
multiply the matrix by itself and stores the result in $k^2$-tree representation in `m2.k2`.
```
k2sparse.x -d m2.k2 -o m2.txt
```
decompress the $k^2$-tree representation of `m2.k2` producing the list of nonzero elements 
in textual form, one nonzero per row. 


### The uncompressed matrix formats

Since all matrices are binary, the uncompressed format consists of a text file containing only the positions of the nonzero entries. Each line should contain the row and column indexes of a single entry written in decimal and separated by a whitespace character.  The same entry should not appear twice, but the order of entries can be arbitrary. Indexes are 0-based see the file `t8.txt` for an 8x8 example. 

For small sizes the tool `asc2sparse.py` can be used to obtain a more human readable representation. For example typing
```
asc2sparse.py -r -o t8.asc t8.txt 
```
produces the file `t8.asc` containing the dense ascii representation of the matric:
```
00000111
00000010
00000011
00000101
11110010
11110001
11110001
11110010
```


### Compression to/from $k^2$-tree format 

The executable `k2sparse.x` is used to compress and decode a single textual matrix to/from a k2tree representation; run it without arguments to get basic usage instructions
```
Usage:
	  k2sparse.x [options] infile

Tool to (de)compress boolean matrices in sparse text format (one line per entry)
to k2 compressed Plain Depth First format

Options:
	-d      decompress to text one line per entry format
	-n      do not write the output file, only show stats
	-o out  outfile name (def. compr: infile.k2, decompr: infile.txt)
	-s S    matrix actual size (def. largest index+1) [compression only]
	-m M    minimatrix size (def. 2) [compression only]
	-i info infile subtree info file [decompression only] (def. None)
	-I info infile backpointers file [decompression only] (def. None)
	-r size rank block size for backpointres (def. 64)
	-x      do not compact all 1's submatrices [compression only]
	-c      compress->decompress->check
	-h      show this help message
	-v      verbose

Default action is to compress filename to filename.k2
```

Invoked without the `-d` or `-c` switch `k2sparse.x` compresses a textual matrix to the k2 format. 

When invoked with `-d` the input file must be a k2 matrix which is then expanded into the textual format. 

When invoked with `-c` the program compresses the input matrix, then decompresses it and verify that the decompressed matrix matches the original matrix. Since the order of the entries in the textual file is arbitrary the verification involves sorting and searching and is done invoking the tool `matrixcmp.x`.

The option `-m` can be used only in compression, currently only with one of the two values `2` or `4`. This parameter is the size of the matrices stored at the leaf of the $k^2$-tree (except the leaves representing the submatrices of all 1's which can be of any size). A large leaf size usually yields larger files but improves the running time for the aritmetic operations over the matrices (as the tree is shallower). The value of the parameter `-m` is stored in the k2 compressed file so it does not have to be provided for decompression.
Other command line options will be explained after we discuss the different compression formats. 

By default the input matrix is assumed to be of size 1+(largest index in the input file). The option `-s` can be used to force the size of the input matrix to a specific (larger) value. Because of the algorithm used to compress textual matrices, currently the largest admissible matrix size is $2^{32}$; this limitation can be removed if needed using a slightly more complex compression algorithm. 





## Product of matrices in k2 format

The executable `k2mult.x` can be used to multiply two compressed matrices in k2 format. The matrices must have the same size and must have been compressed with the same `-m` parameter. The output is still in k2 format.
```
Usage:
	  k2mult.x [options] infile1 infile2

Options:
	-n        do not write output file, only show stats
	-o out    outfile name (def. infile1.prod)
	-i info   infile1 subtree info file
	-j info   infile2 subtree info file
	-I info   infile1 backpointers file
	-J info   infile2 backpointers file
	-t size   rank block size for k2 compression (def. 64)
	-e        compute subtree info on the fly (def. no)
	-x        do not compact new 1's submatrices in the result matrix
	-q        use a single copy when squaring a matrix
	-c        check multiplication (O(n^3) time and O(n^2) space!)
	-h        show this help message
	-v        verbose

Multiply two compressed matrices stored in infile1 and infile2
```
When invoked with `-c`, after computing the product in k2 format, the program uncompresses the input matrices and the product and verify that the uncompressed product is identical to the product computed with the traditional $O(n^3)$ time algorithm applied to the uncompressed inputs. For large matrices this verification can be slow and space consuming. 

For an explanation of the `-e`, `-i` and `-j` options see section *Enriched format* below. 


### Example:

The following sequence of compression, matrix multiplication, decompression and display operations. 
```
k2sparse.x t8.txt
k2mult.x t8.txt.k2 t8.txt.k2 -o t8sq.k2
k2sparse.x -d t8sq.k2 -o t8sq.txt
asc2sparse.py -r -o t8sq.asc t8sq.txt
cat t8sq.asc
```
should eventually display the input matrix squared:
```
11110011
11110001
11110011
11110011
11110111
11110111
11110111
11110111
```

## Compression formats


### Enriched format

In order to speedup operations on compressed matrices in k2 format it is possible to use some extra information on its largest subtrees. This information must be computed with `k2subtinfo.x`, for example:
```bash
k2subtinfo.x -vc m.k2 -o m.k2.sinfo 
```
computes the information for the compressed matrix `m.k2` and stores it in `m.k2.sinfo`. By default the information is computed for the subtrees whose size is at least the square root of the size of the k2 tree representing the whole matrix (here size is measured in number of nodes). 
Use the option `-N limit` to compute the info only for subtrees of size larger than `limit`. Alternatively use the option `-M` set the limit of to a fraction of the square root of the number of nodes, for example for a tree with 1.000.000 nodes, using `-M 0.2` will set the limit to 200 nodes. As an alternative to `-N` and `-M` one can use the option `-D d` to compute the info only for the subtree at depth up to `d`. 

At the moment this additional information can be used only to speed-up matrix-matrix multiplication. For example write
```bash
k2mult.x m1.k2 m2.k2 -i m1.k2.sinfo -j m2.k2.sinfo 
```
to compute the product `m1.k2` times `m2.k2` using the extra information `m1.k2.sinfo` for the compressed matrix `m1.k2`and `m2.k2.sinfo` for the compressed matrix `m2.k2` (the information can be specified also for only one of the two input matrices). 

The option `-e` enables the "on the fly" computation of the subtree information: when the multiplication algorithm reaches a subtree for which no subtree information is available, the information is computed and later discarded.  

NOTE: to see a real speed improvement in the matrix multiplication algorithm using subtree size info at the moment it is necessary to use the release version (`make release`) since the default version still performs many redundant checks. 



### Compressed k2-tree

The executable `k2cpdf.x` is used to compress a $k^2$-tree representation using subtree compression; run it without arguments to get basic usage instructions
```
Usage:
	  k2cpdf.x [options] infile

Options:
	-b      block size for rank 0000 operation (def. 64)
	-t      smallest subtree size to be removed in bits (def. 32)
	-c      check number of ones in the compressed matrix
	-n      do not write the output file, only show stats
	-h      show this help message
	-v      verbose

Compress a k2 tree exploiting the presence of identical subtrees.
Compute and store in separates files the compressed tree and
its auxiliary information (pointers information)
```

When invoked will generate three new files in the same directory as `infile`:

* `infile.ck2`: file containing the compressed k2-tree: repeated subtree are  represented by the special node `0000`.
* `infile.ck2.p`: a binary file holding a `uint64_t` array. Each value is the destination of a special pointer node `0000`.  


When invoked with `-c` will check that the compressed representation contains the same amount of non-zero elements as the original k2tree representation.

When invoked with `-n` will not write any file, only compress and show statistic.


### Enriched compressed format

To compute subtree information for compressed k2-tree use `k2subtinfo.x` with the `-p` options as follows
```bash
k2subtinfo.x -p -vc m.ck2
```
This command creates a file `m.ck2.sinfo` containing the subtree information..

To compute the matrix product using the enriched compressed format:
```bash
k2mult.x -v -i a.ck2.sinfo -I a.ck2.p  -j b.ck2.sinfo -J b.ck2.p  a.ck2  b.ck2
```



## Pagerank computation 

The Pagerank computation is ideal for testing the speed of the matrix-vector product in a real-world scenario.
Assuming that the input (web) matrix is given in mtx format it is first necessary to preprocess 
it using the `mtx2rowm` tool that after the conversion provides some minimal instructions to
compress the input matrix (using `k2sparse.x` or `k2blockc.py`) and later compute the Pagerank
vector (using `k2pagerank.x`) possibly using multiple threads. 



## Matrices represented as bitarrays

The library also contains the code for compressing and operating on boolean matrices using a bitarray, ie using one bit per entry plus a small overhead. To make the conversion between the two compressed formats very simple, the callable functions (whose prototypes are in `k2.h` and `b128.h`) have the same names. Hence, a program using the k2 format can be transformed into one using the bitarray format by redefining a few constants. See the use of the `B128MAT` compilation constant in the source files `k2bbm.c` and `k2mult.c` and in the `makefile`. Creation of bitarray matrices is currently not supported for textual input matrices. Since bitarray representation does not take advantage of sparsity, the largest supported size is $2^{30}$. 

The programs `b128sparse.x`, `b128bbm.x` and `b128mult.x` work exactly like  `k2sparse.x`, `k2bbm.x` and `k2mult.x` except that they use the bitarray representation instead of the k2 format. 



## Additional tools 


The program `k2showinfo.x` display statics on the k2-compressed files passed on the command line.

The script `submatrix.py` can be used to extract a square submatrix form a matrix in textual form.

The files `k2test.sh`, `k2btest.sh`, `k2square.sh` are bash script designed to test `k2sparse.x`, `k2bbm.x` and `k2mult.x` on a set of input files.  


## Tools no longer supported

### Product of bbm matrices with openmp

The program `bbmmult.x` computes the product of two `.bbm` matrices using `openmp` to speedup the computation. This tool has been provided mainly as an alternative to the `-c` option to check the correctness of `k2mult.x` and `b128mult.x`.
