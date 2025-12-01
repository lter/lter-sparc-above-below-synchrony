#!/bin/bash
#SBATCH --job-name=NEON_REPLACEME_ASV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=gzahn@wm.edu
#SBATCH -o /sciclone/home/gzahn/logs/NEON_DADA_REPLACEME_%j.out
#SBATCH -e /sciclone/home/gzahn/logs/NEON_DADA_REPLACEME_%j.err




Rscript /sciclone/scr10/gzahn/Synchrony/jobs/REPLACEME_asv_2020.R
