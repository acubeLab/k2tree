#!/usr/bin/env bash

# test squaring for a set of matrices already in k2 format 
# comparing results with precomputed sha1sum values

# exit immediately on error
set -e

if [ $# -lt 2 ]
then
  echo "Usage:"
  echo "         $0 dir file1 [file2 ...]"
  echo
  echo "Test squaring with subtinfo, backpointers, and dynamic subtrees"
  echo "for file1.k2 file2.k2 ... all from directory dir"
  echo "Also compute sha1sums and compare with those in prod.sha1sum" 
  echo
  echo "Sample usage (better with nohup):"
  echo "         $0 web cnr80k eu80k goo120k"        
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
  echo "====== squaring with dynamic subtrees: $f"
  $timecmd -f"$tf" ./k2mult.x -e -xq $dir/$f.k2  $dir/$f.k2 
  sha1sum --ignore-missing -c $dir/prod.sha1sum

  echo "==== adding subtree information to $f ==="
  ./k2subtinfo.x $dir/$f.k2
  echo "====== squaring with subtree info: $f"
  $timecmd -f"$tf" ./k2mult.x -xq $dir/$f.k2  $dir/$f.k2 -i $dir/$f.k2.sinfo
  sha1sum --ignore-missing -c $dir/prod.sha1sum
  echo "====== squaring with subtree info + dynamic subtrees: $f"
  $timecmd -f"$tf" ./k2mult.x -xq $dir/$f.k2  $dir/$f.k2 -i $dir/$f.k2.sinfo
  sha1sum --ignore-missing -c $dir/prod.sha1sum

  echo "==== compressing subtrees on $f"
  $timecmd -f"$tf" ./k2cpdf.x $dir/$f.k2
  echo "====== squaring with compressed subtrees: $f"
  $timecmd -f"$tf" ./k2mult.x -xq $dir/$f.ck2  $dir/$f.ck2 -I $dir/$f.ck2.p -o $dir/$f.k2.prod
  sha1sum --ignore-missing -c $dir/prod.sha1sum

  echo "==== adding subtree information to $f.ck2"
  $timecmd -f"$tf" ./k2subtinfo.x $dir/$f.ck2
  echo "====== squaring with compressed subtrees + subtree info: $f"
  $timecmd -f"$tf" ./k2mult.x -xq $dir/$f.ck2  $dir/$f.ck2 -I $dir/$f.ck2.p -i $dir/$f.ck2.sinfo -o $dir/$f.k2.prod
  sha1sum --ignore-missing -c $dir/prod.sha1sum
  echo "====== squaring with compressed subtrees + subtree info + dynamic subtrees:$f"
  $timecmd -f"$tf" ./k2mult.x -e -xq $dir/$f.ck2  $dir/$f.ck2 -I $dir/$f.ck2.p -i $dir/$f.ck2.sinfo -o $dir/$f.k2.prod
  sha1sum --ignore-missing -c $dir/prod.sha1sum
  # echo "==== deleting $f.k2.prod"
  # rm -f $dir/$f.k2.prod

done
