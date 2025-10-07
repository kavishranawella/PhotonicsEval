#!/bin/bash

#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH -p cpu
#SBATCH --job-name=frozenk_nL_pde
#SBATCH --mem=250000
#SBATCH --output=output_frozenk_nL_pde.txt

./PHOTONICS/frozenk_nL_pde.out