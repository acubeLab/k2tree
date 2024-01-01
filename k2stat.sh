#!/usr/bin/env bash

# compress decompress and chek all files passed one the command line


if [ $# -le 2 ]
then
  echo "Usage:"
  echo "         $0 compr-options extension file1 [file2 ...]"
  echo
  echo "Report the compression statistics for a single k2comp variant"
  echo "on all input files (wildcards in file names are ok)"
  echo
  echo "Sample usage:"
  echo "         $0 \"-m4\" \".k4\" data/*.bbm"        
  exit
fi

options=$1
ext=$2
shift 2
echo $options
for f in "$@"
do 
  echo ">>>>>>>> File: $f"
  echo "==== compression"
  k2comp.x $options -o $f$ext  $f
  ls -l $f$ext
  echo "==== decompression" 
  k2comp.x -v -d $f$ext
  echo "==== check"
  cmp $f $f$ext.bbm
done
