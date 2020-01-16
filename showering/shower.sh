#!/bin/bash

# creating folders
mkdir -p ggh-csv
mkdir -p vbf-csv
mkdir -p vh-csv

#showering
./myexample.exe signal.cmnd ggh-csv ggh-lhe/*.lhe >> ggh-csv/shower.log
./myexample.exe signal.cmnd vbf-csv vbf-lhe/*.lhe >> vbf-csv/shower.log
./myexample.exe signal.cmnd vh-csv vh-lhe/*.lhe >> vh-csv/shower.log

