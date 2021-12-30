#!/bin/bash -l

#SBATCH --job-name=polyfreqs_par
#SBATCH -o polyfreqs_par.out
#SBATCH -e polyfreqs_par.err
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --mail-user=yshi@uwsp.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --workdir=/scratch-lustre/wlarson/yue/sturgeon/polyfreqs

enable_lmod
module load R/3.6
module load polyfreqs

Rscript ../scripts/polyfreqs_parallel.r
