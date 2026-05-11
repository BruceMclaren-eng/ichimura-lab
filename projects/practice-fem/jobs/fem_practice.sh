#!/bin/bash
#SBATCH --job-name=fem_practice
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --ntasks=1
#SBATCH --time=00:10:00

./fem_practice
