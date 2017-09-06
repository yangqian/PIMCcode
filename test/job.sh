#!/bin/bash  
NPARTICLE=2 #total number of particles
NPARTICLE1=1 #number of particles of species 1
NPARTICLE2=1 #number of particles of species 2
BEADS=4 #total number of beads
WIGGLEBEADS=2 #number of wiggle beads
BETA=0.5 #inverse temperature
OMEGA=1 #trapping frequency
COL=1 #correlation length
DIL=1.7 #dialation move length
DIS=0.55 #single slice move length
SHIFT=0 #disable whole path move (negligible)
RUNPATH=../fermi/run #path of the executable
NPROC=4 #total number of processors
#for actual run, remove the --seed=12345678
rm fpath.save*
echo "Starting the simulation from scratch, only calculating scalar observables."
echo "Random number seed 12345678; should be changed for an actual simulation."
echo "----------------------------------------------------------------"
BLOCKLENGTH=400 #single block length
args="--beta=$BETA --blocklength=$BLOCKLENGTH --Nbeads=$BEADS --Nwiggle=$WIGGLEBEADS --Nparticle1=$NPARTICLE1 --Nparticle2=$NPARTICLE2 --col=$COL -u --initialize=1 -d$DIL --shift=$SHIFT --displacement=$DIS --centerofmassnon"
time mpirun -n $NPROC $RUNPATH $args --energyonly=1 --seed=12345678
echo "----------------------------------------------------------------"
BLOCKLENGTH=1000 #single block length
args="--beta=$BETA --blocklength=$BLOCKLENGTH --Nbeads=$BEADS --Nwiggle=$WIGGLEBEADS --Nparticle1=$NPARTICLE1 --Nparticle2=$NPARTICLE2 --col=$COL -u --initialize=1 -d$DIL --shift=$SHIFT --displacement=$DIS --centerofmassnon"
echo "Starting a longer simulation with saved configuration."
echo "Random number seed 12345678."
echo "----------------------------------------------------------------"
time mpirun -n $NPROC $RUNPATH $args --seed=12345678
echo "----------------------------------------------------------------"
echo "run cat col1.dat|../b4 to get a detailed analysis of the second estimator"
echo "list the mean, error, error contribution by sign, and statistics"
echo "----------------------------------------------------------------"
cat col1.dat|../b4
echo "----------------------------------------------------------------"
echo "scaled single particle distribution stored as fwave1.dat"
echo "header reads x average std 0 1"
echo "means position, histogram average, standard deviation,"
echo "distribution for particle 0, distribution for particle 1"
echo "unscaled single particle distribution stored as fwave1d.dat"
echo "Change corresponding code in calc.c to adjust correlation choice."
echo "Change corresponding code in initialize.c to adjust estimator choice."
