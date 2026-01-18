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

# compute sa
for f in "$@"
do 
  echo "====== executing product with on the fly computation for file $f"
  ./k2mult.x $dir/$f.k2  $dir/$f.k2 
  sha1sum --ignore-missing -c $dir/prod.sha1sum
  echo "==== adding subtree information ==="
  ./k2subtinfo.x $dir/$f.k2
  echo "====== executing product with subtree info on file $f"
  ./k2mult.x -q $dir/$f.k2  $dir/$f.k2 -i $dir/$f.k2.sinfo
  sha1sum --ignore-missing -c $dir/prod.sha1sum
  echo "==== compressing subtrees on file $f"
  ./k2cpdf.x $dir/$f.k2
  echo "====== executing product with compressed subtree on file $f"
  ./k2mult.x -q $dir/$f.ck2  $dir/$f.ck2 -I $dir/$f.ck2.p -o $dir/$f.k2.prod
  sha1sum --ignore-missing -c $dir/prod.sha1sum
  echo "==== adding subtree information for $f.ck2"
  ./k2subtinfo.x $dir/$f.ck2
  echo "====== executing product with compressed subtree and info on file $f"
  ./k2mult.x -q $dir/$f.ck2  $dir/$f.ck2 -I $dir/$f.ck2.p -i $dir/$f.ck2.sinfo -o $dir/$f.k2.prod
  sha1sum --ignore-missing -c $dir/prod.sha1sum
done

# delete prod file 
for f in "$@"
do 
  echo "==== deleting $f.k2.prod"
  rm -f $dir/$f.k2.prod
done
