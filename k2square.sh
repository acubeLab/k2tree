#!/usr/bin/env bash

# compress decompress and chek all files passed one the command line

if [ $# -le 2 ]
then
  echo "Usage:"
  echo "         $0 file1 [file2 ...]"
  echo
  echo "Compute the square of all input k2-mats (wildcards in file names are ok)"
  echo "reporting running time and peak memory usage"
  echo
  echo "Sample usage:"
  echo "         $0 data/*.k2"        
  exit
fi


for f in "$@"
do 
  echo ">>>>>>>> Squaring $f"
  /usr/bin/time -f"Command: %C\nE(secs):%e Mem(kb):%M" k2mult.x $f $f
done
