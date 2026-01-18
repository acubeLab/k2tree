# Compression formats and operations 


## Plain Depth Dirst 

Plain Depth First representation with compression of all 1 submatrices

### compression

```
k2sparse.x m.txt -1 -o m.k2
```

### squaring

```
k2mult.x -q m.k2 -o msquare.k2
```

### multiplication
```
k2mult.x m1.k2 m2.k2 -o m3.k2
```


## Enriched Depth First 

As above but uses subtree size information for the upper levels

### compression

```
./k2sparse.x m.txt -1 -o m.k2
./k2subtinfo.x m.k2 
```

### squaring

```
k2mult.x -q m.k2 -i m.k2.sinfo -o msquare.k2
```

### multiplication
```
k2mult.x m1.k2 m2.k2 -i m1.k2.sinfo -j m2.k2.sinfo -o m3.k2
```

### squaring/multiplication with dynamic subtree info

To compute the subtree info on the fly for small subtrees use option `-e`; for example
```
k2mult.x -e -q m.k2 -i m.k2.sinfo -o msquare.k2
```


## Compressed Depth First 

Compress the $k^2$ tree by detecting identical subtrees. No authomatic compression of submatrices of all 1s.

###  Compression

```
./k2sparse.x m.txt -o m.k2
./k2cpdf.x 
./k2subtinfo.x m.ck2 
```


### squaring

```
k2mult.x -q m.ck2 -i m.ck2.sinfo -I m.ck2.p -o msquare.k2
```

### multiplication
```
k2mult.x m1.k2 m2.k2 -i m1.k2.sinfo -I m1.ck2.p -j m2.k2.sinfo -J m2.ck2.p -o m3.k2
```

### squaring/multiplication with dynamic subtree info

AS above, to compute the subtree info on the fly for small subtrees use option `-e`; for example
```
k2mult.x -e -q m.ck2 -i m.ck2.sinfo -I m.ck2.p -o msquare.k2
```


# Available C functions 