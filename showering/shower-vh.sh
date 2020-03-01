#!/bin/bash

# creating folders
mkdir -p vh-csv

# Run showering for every <seed>.lhe file, then save as vh-csv/<seed>.lhe.csv

for f in vh-lhe/*.lhe;
  do 
  echo $'processing file: '$f 
  seed=$(echo $f|sed 's/[^0-9]*//g') 
  ./myexample.exe hj.cmnd vh-csv/$seed $f >> vh-csv/shower.log    
  done



