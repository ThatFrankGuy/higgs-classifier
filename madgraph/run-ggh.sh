#!/bin/bash

# make a new folder for generated stuff, in case it doesnt exist
mkdir -p ggh-run

# generate a new seed.
seed=$RANDOM

# generate script and log file names
scriptnamehj=$PWD$"/ggh-run/mg-script-ggh-hj-"$seed$".txt"
lognamehj=$PWD$"/ggh-run/mg-script-ggh-hj-"$seed$".log"
scriptnamehjj=$PWD$"/ggh-run/mg-script-ggh-hjj-"$seed$".txt"
lognamehjj=$PWD$"/ggh-run/mg-script-ggh-hjj-"$seed$".log"

#see if the files exist, if so, generate a new set

while [ -f $scriptname ]
do
   seed=$RANDOM
   scriptnamehj=$PWD$"/ggh-run/mg-script-ggh-hj-"$seed$".txt"
   lognamehj=$PWD$"/ggh-run/mg-script-ggh-hj-"$seed$".log"
   scriptnamehjj=$PWD$"/ggh-run/mg-script-ggh-hjj-"$seed$".txt"
   lognamehjj=$PWD$"/ggh-run/mg-script-ggh-hjj-"$seed$".log"
done

# create a script file containing unique seed
echo -e "launch ~/scratch/higgs-classifier/MG5/ggh-hj" >> $scriptnamehj
echo -e $"set iseed "$seed >> $scriptnamehj
echo -e "launch ~/scratch/higgs-classifier/MG5/ggh-hjj" >> $scriptnamehjj
echo -e $"set iseed "$seed >> $scriptnamehjj

# setup modules
module unload PrgEnv-intel
module load PrgEnv-gnu

# launch madgraph processes
echo -e $"Start time: "$date >> $lognamehj

/global/homes/f/frankfu/scratch/software/MG5_aMC_v2_6_7/bin/mg5_aMC $scriptnamehj >> $lognamehj
/global/homes/f/frankfu/scratch/software/MG5_aMC_v2_6_7/bin/mg5_aMC $scriptnamehjj >> $lognamehjj

echo -e $"End time: "$date >> $lognamehj
