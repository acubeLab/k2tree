#!/usr/bin/env bash

# test matrix multiplication speed squaring each input matrix

# exit immediately on error
set -e


if [ $# -lt 1 ]
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
  /usr/bin/time -f"Command: %C\nE(secs):%e Mem(kb):%M" k2mult.x -c $f $f
done
