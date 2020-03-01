#!/bin/bash

# creating folders
mkdir -p ggh-hj-csv

# Run showering for every <seed>.lhe file, then save as vh-csv/<seed>.lhe.csv

for f in ggh-hj-lhe/*.lhe;
  do 
  echo $'processing file: '$f 
  seed=$(echo $f|sed 's/[^0-9]*//g') 
  ./myexample.exe hj.cmnd ggh-hj-csv/$seed $f >> ggh-hj-csv/shower.log    
  done



