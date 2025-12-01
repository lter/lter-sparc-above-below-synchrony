#!/bin/bash
#SBATCH --job-name=MG_Assembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=gzahn@wm.edu               # change this!
#SBATCH -o /sciclone/home/gzahn/logs/NEON_download_%j.out # change this!
#SBATCH -e /sciclone/home/gzahn/logs/NEON_download_%j.err # change this!


Rscript /sciclone/scr10/gzahn/Synchrony/download_neon.R
