#!/usr/bin/env bash

# exit immediately on error
set -e

if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "         $0 file1 [file2 ...]"
  echo
  echo "Test A+0 A+I A+A^t on sparse text files file1 file2 ... "
  echo
  echo "Sample usage (better with nohup):"
  echo "         $0 web/*.txt"        
  exit
fi

# A+0 and A+I are computed by k2unary.x
# A+A^t is computed by k2sum.x
# the results are checked by comparing with those obtained 
# by sparsetr.py operating on the sparse representation in the input file



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

  
  echo "--------------------------------------------------------"
  echo "   Working on matrix $f"
  echo "--------------------------------------------------------"

  # testing k2unary.x ie A+0 and A+I
  echo "===== add diagonal to sparse matrix: $f"  
  $timecmd -f"$tf" util/sparsetr.py -d -o $f.sparse.1 $f

  echo "===== compress and check matrix $f in formats k2 and k4"
  $timecmd -f"$tf" ./k2sparse.x -cv -o $f.k2  $f 
  $timecmd -f"$tf" ./k2sparse.x -c -o $f.k4 -m4 $f 

  echo "==== add 0, add I and compute (A+I)^2  ==="
  $timecmd -f"$tf" ./k2unary.x -o $f.2  $f.k2 
  $timecmd -f"$tf" ./k2unary.x -o $f.4  $f.k4 

  echo "==== check k2+0 matrix  ==="  
  $timecmd -f"$tf" ./matrixcmp.x $f.2.0.txt $f
  echo "==== check k2+I matrix  ==="  
  $timecmd -f"$tf" ./matrixcmp.x $f.2.1.txt $f.sparse.1

  echo "==== check k4+0 matrix  ==="  
  $timecmd -f"$tf" ./matrixcmp.x $f.4.0.txt $f
  echo "==== check k4+I matrix  ==="  
  $timecmd -f"$tf" ./matrixcmp.x $f.4.1.txt $f.sparse.1

  # delete some temp file
  rm -f  $f.4.0.txt $f.4.1.txt $f.4.1sq.txt

  # testing k2sum.x A+A^t
  echo "===== transpose sparse matrix: $f"  
  $timecmd -f"$tf" util/sparsetr.py -t -o sparse.tr $f
  # this is A+A^t = A^t + A
  echo "===== make symmetric sparse matrix: $f"  
  $timecmd -f"$tf" util/sparsetr.py -S -o sparse.sym $f

  echo "===== compress transposed matrix $f in formats k2 and k4"
  $timecmd -f"$tf" ./k2sparse.x -o $f.tr.k2  sparse.tr 
  $timecmd -f"$tf" ./k2sparse.x -o $f.tr.k4 -m4 sparse.tr 
  echo "===== compute k2 + k2^t ======="
  $timecmd -f"$tf" ./k2sum.x -o $f.sum.k2 $f.k2 $f.tr.k2 
  echo "===== check k2 + k2^t ======="
  $timecmd -f"$tf" ./k2sparse.x -d -o sparse.sum $f.sum.k2
  $timecmd -f"$tf" ./matrixcmp.x sparse.sum sparse.sym
  echo "===== compute k4^t + k4 ======="
  $timecmd -f"$tf" ./k2sum.x -o $f.sum.k4 $f.tr.k4  $f.k4
  echo "===== check k4^t + k4 ======="
  $timecmd -f"$tf" ./k2sparse.x -d -o sparse.sum $f.sum.k4
  $timecmd -f"$tf" ./matrixcmp.x sparse.sum sparse.sym

  # delete some temp file
  rm -f  sparse.tr $f.k4 $f.tr.k4 $f.sum.k4 
  
  echo "====== compress matrix $f in format k2 + backpointers"
  # $timecmd -f"$tf" ./k2sparse.x -x -o $f.k2  $f 
  $timecmd -f"$tf" ./k2cpdf.x $f.k2
  if [ -e $f.ck2.p ]; then
    echo "== uncompress and test ck2 matrix"
    $timecmd -f"$tf" ./k2sparse.x -d $f.ck2 -I $f.ck2.p -o sparse.ck2
    $timecmd -f"$tf" ./matrixcmp.x $f sparse.ck2
    echo "== unary operations"
    $timecmd -f"$tf" ./k2unary.x -o $f.2 -I $f.ck2.p $f.ck2 
    echo "== check ck2+0 matrix"  
    $timecmd -f"$tf" ./matrixcmp.x $f.2.0.txt $f
    echo "== check ck2+I matrix"  
    $timecmd -f"$tf" ./matrixcmp.x $f.2.1.txt $f.sparse.1
    echo "== sum k2^t + ck2"
    $timecmd -f"$tf" ./k2sum.x -o $f.sum.k2 -J $f.ck2.p $f.tr.k2 $f.ck2
    echo "== check k2 + k2^t"
    $timecmd -f"$tf" ./k2sparse.x -d -o sparse.sum $f.sum.k2
    $timecmd -f"$tf" ./matrixcmp.x sparse.sum sparse.sym
  else
    echo "===== no backpointers for file $f"
  fi

  echo "==== cleaning ==="
  rm -f $f.2.0.txt $f.2.1.txt $f.2.1sq.txt $f.sparse.1
  rm -f $f.ck2 $f.ck2.p
  rm -f $f.sum.k2 $f.tr.k2 $f.k2
  rm -f sparse.sum sparse.sym sparse.ck2

done

# delete tmp files 
# for f in "$@"
# do 
#   echo "==== deleting $f.2.tr.txt"
#   rm -f $f.2.tr.txt $f.2.tr1.txt
# done
