#!/usr/bin/zsh

d=$1
cwd=$(pwd)
cd $d || exit

for f in $( find ./ -name "*.dcd")
do
  fname=$( basename "$f" | cut -f 1 -d '.')
  pdb=protein.pdb
  echo $fname
  mdconvert -o "$fname".xtc -s 5 -t "$pdb" "$fname".dcd
done
cd $cwd





