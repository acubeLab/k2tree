#!/usr/bin/env bash

# exit immediately on error
set -e

# report compression statistics for a set of files
# do not check correctness of compressor: use k2test.sh for that

if [ $# -le 2 ]
then
  echo "Usage:"
  echo "         $0 compr-options outfile_ext file1 [file2 ...]"
  echo
  echo "Report the compression statistics for a single k2sparse variant"
  echo "on all input files (wildcards in file names are ok)"
  echo
  echo "Sample usage:"
  echo "         $0 \"-m4\" \".k4\" data/*.txt"        
  exit
fi

options=$1
ext=$2
shift 2
for f in "$@"
do 
  echo ">>>>>>>> File: $f"
  k2sparse.x $options -vc -o $f$ext $f
  rm -f $f$ext.check
done
