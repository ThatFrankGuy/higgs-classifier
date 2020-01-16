Joshua Lin, Modified and updated by Frank Fu

To compile:
This code uses fastjet 3.3.3, pythia 8.2.4.4 and fastjet-contrib 1.042
first install all packages using instructions given at 
http://home.thep.lu.se/Pythia/
http://fastjet.fr/quickstart.html
https://fastjet.hepforge.org/contrib/ (note: remember to configure --fastjet-config=FILE)
fastjet-contrib's ConstituentSubstractor doesn't work on Cori, and needed to be disabled in makefile. This modified makefile is stored in fastjet-contrib-makefile.
Before compiling, run `source setup-cori.sh` to setup library path (you need to change all paths to your own installation), `make clean` and `make`. 
The compiled file will be called `./myexample.exe`.

To run: 
./myexample.exe signal.cmnd (LOCATION TO STORE EVENTS) (LHE OUTPUT FROM MG5)

A script, `shower.sh`, is written to perform showering for all .lhe files in `/<process>-lhe` and store outputs in `/<process>-csv`.
