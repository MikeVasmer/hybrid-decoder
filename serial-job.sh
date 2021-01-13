#!/bin/bash
#SBATCH --account=def-raymond
#SBATCH --time=12:00:00
#SBATCH --mem=2G
#SBATCH --job-name=hybrid
#SBATCH --output=data/%x-%j.out
#SBATCH --mail-user=mvasmer@pitp.ca
#SBATCH --mail-type=ALL

pwd
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
date

module load python 
module load scipy-stack

file="input/19_10_20a.csv"
dech=80
diw=159
python monte_carlo.py $file $dech $diw

date
