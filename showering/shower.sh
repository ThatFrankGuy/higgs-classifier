#!/bin/bash
# Before running, please run: 
#source /global/project/projectdirs/atlas/scripts/setupATLAS.sh
#setupATLAS -c slc6+batch
#lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"


# creating folders
mkdir -p ggh-csv
mkdir -p vbf-csv
mkdir -p vh-csv

# Run showering for every <seed>.lhe file, then save as vh-csv/<seed>.lhe.csv
for f in vh-lhe/*.lhe;
  do 
  echo $f
  seed=$f | sed 's/[^0-9]*//g' 
  ./myexample.exe signal.cmnd vh-csv/$seed.csv $f >> vh-csv/shower.log    
  done

for f in vbf-lhe/*.lhe;
  do 
  echo $f
  seed=$f | sed 's/[^0-9]*//g' 
  ./myexample.exe signal.cmnd vbf-csv/$seed.csv $f >> vbf-csv/shower.log    
  done

for f in ggh-lhe/*.lhe;
  do 
  echo $f
  seed=$f | sed 's/[^0-9]*//g' 
  ./myexample.exe signal.cmnd ggh-csv/$seed.csv $f >> ggh-csv/shower.log    
  done



