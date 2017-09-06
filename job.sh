#!/bin/bash  
NPARTICLE=2
NPARTICLE1=1
NPARTICLE2=1
BLOCKLENGTH=2
BEADS=4
WIGGLEBEADS=2
BETA=0.5
RANGE=0.02
OMEGA=1
COL=1
PERMUTETRIAL=0
PERMUTEBEADS=8
DIL=1.7
PIGS=0.32
DIS=0.55
POWER=0.335
RUNPATH=~/b45/fermi/run
BLOCKLENGTH=40
args="--pigs=$PIGS --beta=$BETA --blocklength=$BLOCKLENGTH --Nbeads=$BEADS --Nwiggle=$WIGGLEBEADS --order=4 --Nparticle1=$NPARTICLE1 --Nparticle2=$NPARTICLE2 --permutebeads=$PERMUTEBEADS --col=$COL --permutetrials=$PERMUTETRIAL -u --omega=$OMEGA --initialize=1 -d$DIL --shift=0 --displacement=$DIS --centerofmassnon --power=$POWER"
time mpirun $nprocflag $RUNPATH $args >out.pre
