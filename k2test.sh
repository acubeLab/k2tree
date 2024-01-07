#!/usr/bin/env bash

# compress/decompress/check all files passed on the command line 

if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "         $0 file1 [file2 ...]"
  echo
  echo "Compress, decompress, and check all the input files with all the"
  echo "known variants of k2sparse (wildcards in file names are ok)"
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
  k2sparse.x -m2 -o $f.k2 -vc1 $f
  echo "==== compression with 4x4 leaves"
  k2sparse.x -m4 -o $f.k4 -vc $f
  echo "==== compression with 4x4 leaves no all ones"
  k2sparse.x -m4 -o $f.k4 -vc $f
  echo "=== delete .check files ==="
  rm -f $f.k2.check $f.k4.check
done


