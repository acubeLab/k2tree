#!/usr/bin/env bash

# compress $1.bin decompress it and compare it with the original

k2comp.x -v -m 4 -e .k4  $1 
k2comp.x -v -d $1.k4
cmp $1 $1.k4.d

