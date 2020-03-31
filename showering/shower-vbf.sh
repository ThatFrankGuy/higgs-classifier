#!/bin/bash

# creating folders
mkdir -p vbf-csv

# Run showering for every <seed>.lhe file, then save as vh-csv/<seed>.lhe.csv

for f in vbf-lhe/*.lhe;
  do 
  echo $'processing file: '$f 
  seed=$(echo $f|sed 's/[^0-9]*//g') 
  ./myexample.exe hj.cmnd vbf-csv/$seed $f >> vbf-csv/shower.log    
  done



