#!/bin/bash
#SBATCH --job-name=NEON_QC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=gzahn@wm.edu               # change this!
#SBATCH -o /sciclone/home/gzahn/logs/NEON_QC_%j.out # change this!
#SBATCH -e /sciclone/home/gzahn/logs/NEON_QC_%j.err # change this!


module load miniforge3
eval "$(conda shell.bash hook)"
conda activate cutadapt-env


Rscript /sciclone/scr10/gzahn/Synchrony/qc_neon_data.R
