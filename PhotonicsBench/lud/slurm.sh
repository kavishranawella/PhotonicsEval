#!/bin/bash

#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH -p cpu
#SBATCH --job-name=lud
#SBATCH --mem=250000
#SBATCH --output=output_lud.txt

./PHOTONICS/lud.out