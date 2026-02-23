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


### The uncompressed matrix format

Since all matrices are binary, the uncompressed format consists of a text file containing only the positions of the nonzero entries. Each line should contain the row and column indexes of a single entry written in decimal and separated by a whitespace character.  The same entry should not appear twice, but the order of entries can be arbitrary. Indexes are 0-based see the file `t8.txt` for an 8x8 example. 

For small sizes the tool `util/asc2sparse.py` can be used to obtain a more human readable representation. For example typing
```
util/asc2sparse.py -r -o t8.asc t8.txt 
```
produces the file `t8.asc` containing the dense ascii representation of the matrix:
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



## Compression formats

By default the representation consist of the traditional $k^2$-tree with each internal node represented by 4 bits denoting whether each of the four submatrices contains some nonzeros or not. The leaves of the tree store explicitly 2x2 matrices still using 4 bits.  
In addition, since the nibble `0000` can never appear as node or leaf, we use it to denote that the correponding submatrix contains only `1`s: as a result such nodes are leaves. The speciality of this implementation is that the tree node are stored in depth first order; so we call this representation Plain Depth First (PDF) format. 

Using the option `-x` in `k2sparse.x` one can prevent the use of `0000` to denote a submatrix of all `1`s which is therefore encoded with the usual rules, but there is no good reason to do that since representing compactly all `1`s submatrices saves space and running time for some operations. 

Using the option `-m4` in `k2sparse.x` on can force the use of 4x4 submatrices as the tree leaves. This will slightly degrade compression but usually improve the running time (there is one level less to traverse).

One can look at the details of the compressed structure with the command `k2showinfo.x`. For example typing 
```
k2showinfo.x web/cnr80k.k2
```
shows the following information:
```
Matrix file: web/cnr80k.k2
cnr80k.k2:
 matrix size: 80000, leaf size: 2, k2_internal_size: 131072
 Nonzeros: 820959, Nonzero x row: 10.261988
 Levels: 17, Nodes: 391845, Minimats: 288040, 1's submats: 14816
```
We see that the leaves store 2x2 matrices and there are 288040 of them (the are callled Minimats). 
Also there are `14816` submatrices containing all `1`'s which are represented with a single `0000` nibble.
The total numer of nonzero elements of 
The input is an 80000x80000 matrix containing 820959 nonzero elements. Since $k^2$-trees represent 
matrices whose size is a power of 2, the input matrix is embedded in a matrix of size $2^{17}$ = 131072.



### Enriched Depth First  (EDF) format

In order to speedup operations on matrices in PDF format it is possible to store some extra information on its largest subtrees. This information must be computed with `k2subtinfo.x`, for example the command:
```bash
k2subtinfo.x -v m.k2 -o web/cnr80k.k2.sinfo web/cnr80k.k2 
```
computes the information for the compressed matrix `web/cnr80k.k2` and stores it in `web/cnr80k.k2.sinfo`. By default the information is computed for the subtrees whose size is at least the square root of the size of the k2 tree representing the whole matrix (here size is measured in number of nodes). Use the option `-N limit` to compute the info only for subtrees of size larger than `limit`. Alternatively, use the option `-M` set the limit to a fraction of the square root of the number of nodes. For example for a tree with 1.000.000 nodes, using `-M 0.2` will set the limit to 200 nodes. As an alternative to `-N` and `-M` one can use the option `-D d` to compute the info only for the subtree at depth up to `d`. 

At the moment this additional information is used only to speed-up matrix-matrix multiplication. For example write
```bash
k2mult.x m1.k2 m2.k2 -i m1.k2.sinfo -j m2.k2.sinfo 
```
to compute the product `m1.k2` times `m2.k2` using the extra information `m1.k2.sinfo` for the compressed matrix `m1.k2`and `m2.k2.sinfo` for the compressed matrix `m2.k2` (the information can be specified also for only one of the two input matrices). 

The option `-e` for `k2mult.x` enables the "on the fly" computation of the subtree information: when the multiplication algorithm reaches a subtree for which no subtree information is available, the information is computed on the fly and later discarded.  

NOTE: to see a real speed improvement in the matrix multiplication algorithm using subtree size info at the moment it is necessary to use the release version (`make release`) since the default version still performs many redundant checks. 



### Subtree compressed k2-tree (CDF) format


In some cases, the $k^2$-tree representation contains identical subtrees, corresponding to identical submatrices on the input matrix. 
Taking advantage of the Depth First representation we can remove identical subtrees and replace them with a pointer to a previous occurrence of the same subtree. 

The executable `k2cpdf.x` is used to compress a $k^2$-tree representation using subtree compression; run it without arguments to get basic usage instructions
```
Usage:
	  k2cpdf.x [options] infile.k2

Options:
	-o out    outfile name (def. infile1.ck2)
	-I info   infile backpointers file (optional)
	-i info   infile subtree info file (not used)
	-r size   block size for rank 0000 data structure (def. 64)
	-t thrs   smallest subtree size in bits to be removed (def. 32)
	-c        check number of ones in the compressed matrix
	-n        do not write the output file, only show stats
	-h        show this help message
	-v        verbose

Compress a k2 tree exploiting the presence of identical subtrees.
Compute and store in separates files the compressed tree and
its auxiliary information (pointers information)
```

The input file must be a $k^2$-tree representation of a matrix; the output consists of two files:

* `infile.ck2`: containing the compressed k2-tree: pruned repeated subtree are represented by the special node `0000` which has the role of a *pointer node*
* `infile.ck2.p`: a binary file holding a `uint32_t` array; each value is the destination of a pointer node.  

Since nibble `0000` is now used as a pointer node it cannot be used to represent a submatrix containing all `1`s. If `0000` nodes representin all `1`s submatrices are present in the input `infile.k2` they will be replaced by a full representation of 

### Enriched compressed format

To compute subtree information for compressed k2-tree use `k2subtinfo.x` as follows
```bash
k2subtinfo.x -p -vc m.ck2
```
This command creates a file `m.ck2.sinfo` containing the subtree information..

To compute the matrix product using the enriched compressed format:
```bash
k2mult.x -v -i a.ck2.sinfo -I a.ck2.p  -j b.ck2.sinfo -J b.ck2.p  a.ck2  b.ck2
```

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

For an explanation of the `-e`, `-i`, `-j`, `-I`, `-J` options see belows the sections on different compression options. 


### Example:

The following sequence of compression, matrix multiplication, decompression and display operations. 
```
k2sparse.x t8.txt
k2mult.x t8.txt.k2 t8.txt.k2 -o t8sq.k2
k2sparse.x -d t8sq.k2 -o t8sq.txt
util/asc2sparse.py -r -o t8sq.asc t8sq.txt
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

## Transitive closure

The `k2tclosure.x` tool can be used to compute the transitive closure of the graph associated  
to a given binary matrix. The computation consists of the repeated iteration of the formula
$A \gets (A+I)A$ until a fixed point is reached (this occurs after at most a logarithmic number of iterations).





## Pagerank computation 

The Pagerank computation is ideal for testing the speed of the matrix-vector product in a real-world scenario.
Assuming that the input (web) matrix is given in mtx format it is first necessary to preprocess 
it using the `util/mtx2rowm` tool that after the conversion provides some minimal instructions to
compress the input matrix (using `k2sparse.x` or `util/k2blockc.py`) and later compute the Pagerank
vector (using `k2pagerank.x`) possibly using multiple threads. 

For example starting with the web graph [cnr-2000.mtx](https://networkrepository.com/web-italycnr-2000.php) 
we proceed as follows:
```
util/mtx2rowm cnr-2000.mtx
./k2sparse.x -s 325557 cnr-2000.mtx.rowm
./k2pagerank.x -v cnr-2000.mtx.rowm.k2 cnr-2000.mtx.ccount
```
should produce the following output
```
==== Command line:
 ./k2pagerank.x -v cnr-2000.mtx.rowm.k2 cnr-2000.mtx.ccount
Number of nodes: 325557
Number of dandling nodes: 78056
Number of arcs: 3216152
Stopped after 100 iterations, delta=4.12028e-07
Sum of ranks: 1.000000 (should be 1)
Top 3 ranks:
  60597 0.023615
  60595 0.023615
  285152 0.009883
Top: 60597 60595 285152
Elapsed time: 8 secs
```


## Additional tools 

The program `k2showinfo.x` display statics on the k2-compressed files passed on the command line.

The script `submatrix.py` can be used to extract a square submatrix form a matrix in textual form.

The files `k2test.sh`, `k2btest.sh`, `k2square.sh` in the directory `utils` are bash scripts designed to test `k2sparse.x`, `k2bbm.x` and `k2mult.x` on a set of input files.  

## Matrices represented as bitarrays

The library also contains the code for compressing and operating on boolean matrices using a bitarray, ie using one bit per entry plus a small overhead. To make the conversion between the two compressed formats very simple, the callable functions (whose prototypes are in `k2.h` and `b128.h`) have the same names. Hence, a program using the k2 format can be transformed into one using the bitarray format by redefining a few constants. See the use of the `B128MAT` compilation constant in the source files `k2bbm.c` and `k2mult.c` and in the `makefile`. Creation of bitarray matrices is currently not supported for textual input matrices. Since bitarray representation does not take advantage of sparsity, the largest supported size is $2^{30}$. 

The programs `b128sparse.x`, `b128bbm.x` and `b128mult.x` work exactly like  `k2sparse.x`, `k2bbm.x` and `k2mult.x` except that they use the bitarray representation instead of the k2 format. 




## Tools no longer supported

### Product of bbm matrices with openmp

The program `bbmmult.x` computes the product of two `.bbm` matrices using `openmp` to speedup the computation. This tool has been provided mainly as an alternative to the `-c` option to check the correctness of `k2mult.x` and `b128mult.x`.
