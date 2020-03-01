#!/bin/bash
lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"


# creating folders
mkdir -p ggh-hj-csv

# Run showering for every <seed>.lhe file, then save as vh-csv/<seed>.lhe.csv

for f in ggh-hj-lhe/*.lhe;
  do 
  echo $f
  seed=$f | sed 's/[^0-9]*//g' 
  ./myexample.exe hj.cmnd ggh-hj-csv/$seed.csv $f >> ggh-hj-csv/shower.log    
  done



