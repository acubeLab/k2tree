# k2tree binary matrix representation

## Prerequisites 

A modern `gcc` and `make`


## Installation 

Clone/download then `make release`



## Getting started


The executable `k2comp.x` is used to compression and decode a single file. 
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

When invoked with `-c` the program compresses the input matrix, then decompresses it and veryfile that the decompressed amtrix matches the original matrix. 


The option `-m` can be used only in compression, currently only with one of the two values `2` and `4`. This parameter is the size of the matrices stored at the leaf of the k2 tree (except the leaves representing the submatrices of all 1's which can be of any size). A large leaf size usually yields larger files but improves the running time for the aritmetic operations over the matrices (as the tree is shallower). THe valu of the parameter `-m` is stored in the k2 format so it does not have to be provided for decompression.


The executable `k2mult.c` is used to multiply two compressed matrices in k2 format. The matrices must have the same size and must have been compressed with the same `-m` parameter. 

```
Usage:
	  k2mult.x [options] iname1 iname2

Options:
	-n      do not write output file, only show stats
	-e ext  extension for the output file (def. .prod)
	-c      check multiplication (very slow for large matrices!)
	-v      verbose

```

