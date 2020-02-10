#!/bin/bash
# this script extracts all .lhe files from ./<process>/Events/run_##/events.lhe.hz

dir=$PWD

mkdir -p ggh-lhe
mkdir -p vbf-lhe
mkdir -p vh-lhe

# Execute for each folder
for f in vh/Events/*;
  do 
    # Find the line with "iseed", remove comment, then remove all non-numeric, leaving the random seed in place.
    seed=$(grep -hr "iseed" $f/run_*.txt | cut -f1 -d"!" | sed 's/[^0-9]*//g' ) 
    # Unzip .gz file
    gzip -d $f/events.lhe.gz
    # copy out .lhe file and rename to random seed
    cp $f/events.lhe $dir/vh-lhe/$seed.lhe
  done

for f in vbf/Events/*;
  do 
    #[ -d $f ] &&
    seed=$(grep -hr "iseed" $f/run_*.txt | cut -f1 -d"!" | sed 's/[^0-9]*//g' ) 
    ls $f

    gzip -d $f/events.lhe.gz
    cp $f/events.lhe $dir/vbf-lhe/$seed.lhe
  done

for f in ggh-hj/Events/*;
  do 
    #[ -d $f ] &&
    seed=$(grep -hr "iseed" $f/run_*.txt | cut -f1 -d"!" | sed 's/[^0-9]*//g' ) 
    ls $f

    gzip -d $f/events.lhe.gz
    cp $f/events.lhe $dir/ggh-lhe/$seed.hj.lhe
  done

for f in ggh-hjj/Events/*;
  do 
    #[ -d $f ] &&
    seed=$(grep -hr "iseed" $f/run_*.txt | cut -f1 -d"!" | sed 's/[^0-9]*//g' ) 
    ls $f

    gzip -d $f/events.lhe.gz
    cp $f/events.lhe $dir/ggh-lhe/$seed.hjj.lhe
  done