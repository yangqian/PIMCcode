#!/bin/bash  
NPARTICLE=3 #total number of particles
NPARTICLE1=2 #number of particles of species 1
NPARTICLE2=1 #number of particles of species 2
BEADS=8 #total number of beads
WIGGLEBEADS=2 #number of wiggle beads
BETA=0.5 #inverse temperature
COL=48 #correlation length
DIL=1.42 #dialation move length
DIS=0.49 #single slice move length
SHIFT=0.39 #whole path move length
RUNPATH=../../fermi/run #path of the executable
NPROC=6 #total number of processors

BLOCKLENGTH=100 #single block length
ENERGYONLY=1 #skip correlator estimators
args="--beta=$BETA --blocklength=$BLOCKLENGTH --Nbeads=$BEADS --Nwiggle=$WIGGLEBEADS --Nparticle1=$NPARTICLE1 --Nparticle2=$NPARTICLE2 --col=$COL -u --initialize=1 -d$DIL --shift=$SHIFT --displacement=$DIS --centerofmassnon -e=ENERGYONLY"
time mpirun -n $NPROC $RUNPATH $args

ENERGYONLY=0
BLOCKLENGTH=1000 #single block length
args="--beta=$BETA --blocklength=$BLOCKLENGTH --Nbeads=$BEADS --Nwiggle=$WIGGLEBEADS --Nparticle1=$NPARTICLE1 --Nparticle2=$NPARTICLE2 --col=$COL -u --initialize=1 -d$DIL --shift=$SHIFT --displacement=$DIS -e=ENERGYONLY"
time mpirun -n $NPROC $RUNPATH $args

#b2 \approx 0.88893
