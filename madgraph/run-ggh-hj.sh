#!/bin/bash

# make a new folder for generated stuff, in case it doesnt exist
mkdir -p ggh-hj-run

# setup path
MG5path=/global/homes/f/frankfu/scratch/software/MG5_aMC_v2_6_7/bin/mg5_aMC
# generate a new seed.
seed=$RANDOM

# generate script and log file names
scriptname=$PWD$"/ggh-hj-run/mg-script-ggh-hj-"$seed$".txt"
logname=$PWD$"/ggh-hj-run/mg-script-ggh-hj-"$seed$".log"

#see if the files exist, if so, generate a new set

while [ -f $scriptname ]
do
   seed=$RANDOM
   scriptname=$PWD$"/ggh-hj-run/mg-script-ggh-hj-"$seed$".txt"
   logname=$PWD$"/ggh-hj-run/mg-script-ggh-hj-"$seed$".log"
done

# create a script file containing unique seed
echo -e "launch ~/scratch/higgs-classifier/MG5/ggh-hj" >> $scriptname
echo -e $"set iseed "$seed >> $scriptname

# setup modules
module unload PrgEnv-intel
module load PrgEnv-gnu

# launch madgraph processes
echo -e $"Start time: "$date >> $logname

$MG5path $scriptname >> $logname

echo -e $"End time: "$date >> $logname
