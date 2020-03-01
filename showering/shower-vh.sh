#!/bin/bash
lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"


# creating folders
mkdir -p vh-csv

# Run showering for every <seed>.lhe file, then save as vh-csv/<seed>.lhe.csv
for f in vh-lhe/*.lhe;
  do 
  echo $f
  seed=$f | sed 's/[^0-9]*//g' 
  ./myexample.exe vh.cmnd vh-csv/$seed.csv $f >> vh-csv/shower.log    
  done

