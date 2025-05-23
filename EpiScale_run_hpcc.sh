#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --output=compileTest# This affects the print out of the "std::cout" in the script, make sure this is changed for different jobs.
#SBATCH --mail-user="enowa001@ucr.edu"
#SBATCH --mail-type=ALL
#SBATCH --job-name="compileTest"
#SBATCH -p gpu # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load cuda

./bin/runDiscSimulation_M -slurm N02_0
