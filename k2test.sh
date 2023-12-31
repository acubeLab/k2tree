#!/usr/bin/env bash

# compress $1.bbm decompress it and compare it with the original

if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "         $0 file1 [file2 ...]"
  echo
  echo "Compress, decompress, and check all the input files with all the"
  echo "known variants of k2-comp (wildcards in file names are ok)"
  echo
  echo "Sample usage:"
  echo "         $0 data/*.bbm"        
  exit
fi

for f in "$@"
do 
  echo ">>>>>>>> File: $f"
  echo "==== compression with 2x2 leaves"
  k2comp.x -m2 -o $f.k2 -cn $f
  echo "==== compression with 4x4 leaves"
  k2comp.x -m4 -o $f.k4 -cn $f
done


