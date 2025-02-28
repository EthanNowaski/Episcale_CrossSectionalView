#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --output=Test# This affects the print out of the "std::cout" in the script, make sure this is changed for different jobs.
#SBATCH --mail-user=enowa001@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Test"
#SBATCH -p gpu # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load singularity
export SINGULARITY_NV=1
module load centos
module load cuda/7.0

centos.sh "module load cuda/7.0; ./bin/runDiscSimulation_M -slurm N02_0"
