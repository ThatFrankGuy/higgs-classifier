#!/bin/bash
alias python=python2
# make a new folder for generated stuff, in case it doesnt exist
mkdir -p vh-run

# generate a new seed.
seed=$RANDOM

# generate script file names
scriptname=$PWD$"/vh-run/mg-script-vh-"$seed$".txt"

#see if the files exist, if so, generate a new set

while [ -f $scriptname ]
do
   seed=$RANDOM
   scriptname=$PWD$"/vh-run/mg-script-vh-"$seed$".txt"
done

# create a script file containing unique seed
echo -e "launch "$PWD$"/vh" >> $scriptname
echo -e $"set iseed "$seed >> $scriptname
echo -e $"set pt_min_pdg {25:450}" >> $scriptname
echo -e $"set pt_max_pdg {25:1700}" >> $scriptname
echo -e $"set nevents 10" >> $scriptname
echo -e $"Launching VH with seed="$seed

# launch madgraph processes
python2 /home/ffu/softwares/MG5_aMC_v2_6_7/bin/mg5_aMC $scriptname


