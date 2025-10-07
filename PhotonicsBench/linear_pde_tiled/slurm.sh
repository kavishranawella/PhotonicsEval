#!/bin/bash

#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH -p cpu
#SBATCH --job-name=linear_pde_tiled
#SBATCH --mem=250000
#SBATCH --output=output_linear_pde_tiled.txt

./PHOTONICS/linear_pde_tiled.out