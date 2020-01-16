#!/bin/bash

setup_PYTHIA() {
    export PYTHIA8LOCATION=/global/homes/f/frankfu/software/pythia8226
    export LD_LIBRARY_PATH=/global/homes/f/frankfu/software/pythia8226/lib/:$LD_LIBRARY_PATH
    export PYTHIA8=/global/homes/f/frankfu/software/pythia8244
    #export PYTHIA8DATA=${PYTHIA8LOCATION}xmldoc/
}

setup_ROOT() {
    lsetup root
}

setup_fastjet() {
    export FASTJETLOCATION=/global/homes/f/frankfu/software/fastjet-install
    export LD_LIBRARY_PATH=/global/homes/f/frankfu/software/fastjet-install/lib/:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTINCDIR=/global/projecta/projectdirs/atlas/bnachman/code/include/
    export BOOSTLIBLOCATION=/global/projecta/projectdirs/atlas/bnachman/code/lib/
    export LD_LIBRARY_PATH=${BOOSTLIBLOCATION}:$LD_LIBRARY_PATH
}

#setup_ROOT
#setup_boost

setup_PYTHIA
setup_fastjet


