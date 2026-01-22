#!/usr/bin/env bash

# exit immediately on error
set -e

if [ $# -lt 2 ]
then
  echo "Usage:"
  echo "         $0 dir file1 [file2 ...]"
  echo
  echo "Test squaring for file1 file2 ... from dir"
  echo  "Also compute sha1sums and compare with those in prod.sha1sum" 
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

  echo "====== transpose and add diagonal sparse matrix: $f"  
  $timecmd -f"$tf" ./sparsetr.py -t -o sparse.tr $dir/$f
  $timecmd -f"$tf" ./sparsetr.py -td -o sparse.tr1 $dir/$f

  echo "====== compress matrix: $f"
  $timecmd -f"$tf" ./k2sparse.x -1 -o $dir/$f.k2  $dir/$f 
  $timecmd -f"$tf" ./k2sparse.x -1 -o $dir/$f.k4 -m4 $dir/$f 

  echo "==== perform transpose and add diagonal ==="
  $timecmd -f"$tf" ./k2unary.x -o $dir/$f.2  $dir/$f.k2 
  $timecmd -f"$tf" ./k2unary.x -o $dir/$f.4  $dir/$f.k4 

  echo "==== check k2 matrices ==="  
  $timecmd -f"$tf" ./matrixcmp.x $dir/$f.2.tr.txt sparse.tr
  $timecmd -f"$tf" ./matrixcmp.x $dir/$f.2.tr1.txt sparse.tr1

  echo "==== check k4 matrices ==="  
  $timecmd -f"$tf" ./matrixcmp.x $dir/$f.4.tr.txt sparse.tr
  $timecmd -f"$tf" ./matrixcmp.x $dir/$f.4.tr1.txt sparse.tr1


done

# delete tmp files 
# for f in "$@"
# do 
#   echo "==== deleting $f.2.tr.txt"
#   rm -f $dir/$f.2.tr.txt $dir/$f.2.tr1.txt
# done
