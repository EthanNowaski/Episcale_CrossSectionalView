#!/bin/bash  

#$ -N "TestingRun"       # Specify job name
#$ -M enowaski@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -q gpu             # Specify queue
#$ -l gpu_card=1         
#$ -pe smp 1        # Specify parallel environment and legal core size

WORK_DIR=$(pwd)
cd $WORK_DIR 

module load cuda

$WORK_DIR/bin/runDiscSimulation_M -slurm N02_0
