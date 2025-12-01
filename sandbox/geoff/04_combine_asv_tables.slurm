#!/bin/bash
#SBATCH --job-name=NEON_COMBINE_ASVs
#SBATCH --nodes=1
#SBATCH -w hi03
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=600
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=gzahn@wm.edu
#SBATCH -o /sciclone/home/gzahn/logs/NEON_COMBINE_%j.out
#SBATCH -e /sciclone/home/gzahn/logs/NEON_COMBINE_%j.err


Rscript /sciclone/scr10/gzahn/Synchrony/combine_asv_tables.R



