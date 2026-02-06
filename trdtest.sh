#!/usr/bin/env bash

# exit immediately on error
set -e

if [ $# -lt 2 ]
then
  echo "Usage:"
  echo "         $0 dir file1 [file2 ...]"
  echo
  echo "Test unary operations on sparse text files file1 file2 ... from dir"
  echo
  echo "Sample usage (better with nohup):"
  echo "         $0 web cnr80k"        
  exit
fi

# get directory
dir=$1
shift 1

echo "Experiment run by $(whoami) at $(date) on $(hostname)"

# compile if necessary 
make

# timing command
# timecmd='/usr/bin/time -f"Command: %C\nS:%S U:%U E:%e Mem(kb):%M\n"'
timecmd=/usr/bin/time
tf="Command: %C\nS:%S U:%U E:%e Mem(kb):%M\n" 

# compute sa
for f in "$@"
do 

  # note tests on matrix transpose removed!

  echo "====== add diagonal to sparse matrix: $f"  
  $timecmd -f"$tf" ./sparsetr.py -d -o sparse.1 $dir/$f

  echo "====== compress matrix: $f in formats k2 and k4"
  $timecmd -f"$tf" ./k2sparse.x -1 -o $dir/$f.k2  $dir/$f 
  $timecmd -f"$tf" ./k2sparse.x -1 -o $dir/$f.k4 -m4 $dir/$f 

  echo "==== add diagonal and compute (A+I)^2==="
  $timecmd -f"$tf" ./k2unary.x -o $dir/$f.2  $dir/$f.k2 
  $timecmd -f"$tf" ./k2unary.x -o $dir/$f.4  $dir/$f.k4 

  echo "==== check k2 matrix ==="  
  $timecmd -f"$tf" ./matrixcmp.x $dir/$f.2.1.txt sparse.1

  echo "==== check k4 matrix ==="  
  $timecmd -f"$tf" ./matrixcmp.x $dir/$f.4.1.txt sparse.1


done

# delete tmp files 
# for f in "$@"
# do 
#   echo "==== deleting $f.2.tr.txt"
#   rm -f $dir/$f.2.tr.txt $dir/$f.2.tr1.txt
# done
