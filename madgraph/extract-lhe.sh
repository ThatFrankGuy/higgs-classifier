#!/bin/bash
# this script extracts all .lhe files from ./<process>/Events/run_##/events.lhe.hz

mkdir -p ggh-lhe
mkdir -p vbf-lhe
mkdir -p vh-lhe

cp ggh-hj/Events/*/*.lhe.gz ggh-lhe
cp ggh-hjj/Events/*/*.lhe.gz ggh-lhe
cp vbf/Events/*/*.lhe.gz vbf-lhe
cp vh/Events/*/*.lhe.gz vh-lhe

tar -xvf ggh-lhe/*.gz
tar -xvf vbf-lhe/*.gz
tar -xvf vh-lhe/*.gz

#rm ggh-lhe/*.gz
#rm vbf-lhe/*.gz
#rm vh-lhe/*.gz
