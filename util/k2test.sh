#!/usr/bin/env bash

# compress/decompress/check all files passed on the command line 

# exit immediately on error
set -e


if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "         $0 file1 [file2 ...]"
  echo
  echo "Compress, decompress, and check all the input files with the"
  echo "plain variants of k2sparse (not including subtree compression"
  echo "with backpointers). Wildcards in file names are ok"
  echo
  echo "Sample usage:"
  echo "         $0 data/*.txt"        
  exit
fi

for f in "$@"
do 
  echo ">>>>>>>> File: $f"
  echo "==== compression with 2x2 leaves"
  k2sparse.x -m2 -o $f.k2 -vc $f
  echo "==== compression with 2x2 leaves no all_ones"
  k2sparse.x -m2 -o $f.k2 -vcx $f
  echo "==== compression with 4x4 leaves"
  k2sparse.x -m4 -o $f.k4 -vc $f
  echo "==== compression with 4x4 leaves no all ones"
  k2sparse.x -m4 -o $f.k4 -vcx $f
  echo "=== delete .k2 .k4 and .check files ==="
  rm -f $f.k2.check $f.k4.check $f.k2 $f.k4
done


