# Introduction

This file contains instructions to perform Monte Carlo simulations for three main types of higgs production events.

Samples are generated for the three following higgs production processes:

1. Gluon-gluon fusion (ggh) (the two subprocesses aren't compatible and need to be run separately)
```
import model loop_sm
generate p p > H j [QCD]
output ggh-hj
```
```
import model loop_sm
generate p p > H j j QED<=2 [QCD]
output ggh-hjj
```
        
2. Vector boson fusion (VBF)
```
import model loop_sm
generate p p > H j j QED=3 [QCD]
output vbf
```

3. Associated production of a Higgs boson and a W or Z boson (VH)
```
import model loop_sm
define v = w+ w- z
generate p p > H v [QCD] @0
add process p p > H v j [QCD] @1
output vh
```

The folder `cards/` contains all the necessary cards to generate 10000 candidate events with 200GeV < pT_higgs < 2000GeV. (since the regime is not resolved, the resulting jet should consist of both bb~ jets and the entire pT of Higgs)

The event selection will be done during parton showing.


`run-<process>.sh` produces 10K sample for a higgs process using a new random seed. All `run-<process>.sh` shares the same MG5 executables, and therefore cannot be exucuted simutaneously by multiple slurm nodes. The log and commands files will be saved in `<process>-run/`. Please change the madgraph program path before using.

`run-<process>.sh` produces 10K sample for a higgs process using a new random seed. It uses temprary madgraph copy, and can be executed simutaneously by multiple slurm nodes. The log and madgraph commands will be saved in `<process>-run/`. Please change the madgraph path before using.

`cori-job-all.sh` runs each madgraph script 50 times to generate 500k samples for each process. For use on Cori Slurm submission.

`extract-lhe.sh` copies all event-level .lhe files from `<process>/` into `<process>-lhe/` and renames them to `<random seed>.lhe` 

**IN ORDER FOR ALL .sh SCRIPTS TO WORK, THE 4 PROCESSES' DIRECTORIES SHOULD BE COPIED INTO A SAME PARENT FOLDER AND RENAMED vbf/, vh/, ggh-hj/ and ggh-hjj/ RESPECTIVELY. ALL .sh SCRIPTS SHOULD BE COPIED TO THE PARENT FOLDER BEFORE EXECUTING.**

Below are the instructions to manually produce higgs classification samples, and to produce samples faster using Cori jobs. 

# Manual instruction

## 0. Installing FastJet
Make sure you have FastJet installed, or things will go wrong at step two. Instruction is [here](http://fastjet.fr/quickstart.html). Once you finish, go into MG5 and execute `set fastjet <installation_path>` to link them.


## 1. Generate Feynman diagrams
Execute codes for each process, then do `output <folder_name>` in MG5 (`<MG5 path>/bin/mg5_aMC`). 
MG5 will iterate possible subprocesses and generate a folder containing their cross sections and Feynman diagrams. 


## 2. Generate events 
Execute `launch <folder_name>` in MG5. It'll ask you about the setting to run events in, like this:
```
The following switches determine which programs are run:
/================== Description ===================|================= values =================|========== other options ===========/
| 1. Type of perturbative computation              |            order = NLO                   |     LO                             |
| 2. No MC@[N]LO matching / event generation       |      fixed_order = OFF                   |     ON                             |
| 3. Shower the generated events                   |           shower = OFF                   |     PYTHIA6Q|PYTHIA6PT|HERWIG6     |
| 4. Decay onshell particles                       |          madspin = OFF                   |     ON|onshell                     |
| 5. Add weights to events for new hypp.           |         reweight = Not Avail.            |     Please install module          |
| 6. Run MadAnalysis5 on the events generated      |      madanalysis = Not Avail.            |     Please install module          |
/==================================================================================================================================/
```   
Or this:
```
/=================== Description ===================|============= values ==============|======== other options ========/
| 1. Choose the shower/hadronization program        |        shower = Not Avail.        |     Please install module     |
| 2. Choose the detector simulation program         |      detector = Not Avail.        |     Please install module     |
| 3. Choose an analysis package (plot/convert)      |      analysis = Not Avail.        |     Please install module     |
| 4. Decay onshell particles                        |       madspin = OFF               |     ON|onshell                |
| 5. Add weights to events for new hypp.            |      reweight = Not Avail.        |     Please install module     |
/=======================================================================================================================/
```
I don't know the reason why they are different for ggh and VBF (Asked Ben and he doesn't know too) but it shouldn't matter. 
Disable shower, madspin, detector, reweight and analysis. Set order to NLO and fixed order to off.

The it will ask you to edit cards. I only edited /Cards/run_card.dat to change event number and add pT cut to Higgs.
Since the run_card.dat will be different for all processes, you'll need to make the edit manually:
    Change `10000 = nevents` to any number you want 
    Change

``` 
  {} = pt_min_pdg ! Min pT for a massive particle
  {} = pt_max_pdg ! Max pT for a massive particle
```

    to

```
  {25:200} = pt_min_pdg ! Max pT for a massive particle
  {25:2000} = pt_max_pdg ! Max pT for a massive particle
```


Some VBF cross sections cannot pass pole check, so for VBF, edit `/Cards/FKS_params.dat`, 
set `#IRPoleCheckThreshold` and `#PrecisionVirtualAtRunTime` to `-1d0`.

Everything should run fine from that. The .lhe files should be in `/Events/run_##`. Note that every time you launch the folder with
different settings MG5 creates a new run_## in events instead of overwritting everything, so take note of the run number you are at.
Or just delete everything in the folder and MG5 will start at run_01 again.

# Cori instruction

## 1. Generate Feynman diagrams
Execute codes for each process, then do `output <folder_name>` in MG5 (`/MG5_aMC_v2_6_6/bin/mg5_aMC`). 
MG5 will iterate possible subprocesses and generate a folder containing their cross sections and Feynman diagrams. 

## 2. Prepare cards
Copy cards in each subfolder of `cards/` into the `<process folder>/Cards/` folder of each corresponding output.

## 3. Run jobs (Cori)
Scripts are provided for the submission, preparation and execution of each process' simulation.

`run-<process>.sh` starts a simulation on the current Cori node using a common MG5 copy. It cannot be run in parallel.
`run-<process>-MG5.sh`  starts a simulation on the current Cori node using a fresh temporary MG5 copy for each simulation. It can be run in parallel. 
`cori-job-<process>.sh` starts a Cori Haswell job.

The scripts should be put in the same folder as MG5 outputs.
Please change the madgraph path when you execute the scripts.

# What's next?

Once the events are generated, go find and copy out the .lhe file from your latest run using `extract-lhe.sh`. Then go use stuff inside `../shower/`.
