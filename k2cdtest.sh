#!/usr/bin/env bash

# compress $1.bin decompress it and compare it with the original

k2comp.x $1 -v
k2comp.x -d $1 -i .k2x -o .k2d -v
cmp $1.bin $1.k2d

