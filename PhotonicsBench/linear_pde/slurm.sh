#!/bin/bash

#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH -p cpu
#SBATCH --job-name=linear_pde
#SBATCH --mem=250000
#SBATCH --output=output_linear_pde.txt

./PHOTONICS/linear_pde.out