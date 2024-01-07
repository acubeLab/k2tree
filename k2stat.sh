#!/usr/bin/env bash

# compress decompress and chek all files passed on the command line
# Obsolete since the introduction of the -c option

if [ $# -le 2 ]
then
  echo "Usage:"
  echo "         $0 compr-options outfile_ext file1 [file2 ...]"
  echo
  echo "Report the compression statistics for a single k2sparse variant"
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
  k2sparse.x $options -v -o $f$ext $f
  rm -f $f$ext.check
done
