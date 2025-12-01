#!/bin/bash
#SBATCH --job-name=NEON_SINTAX
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=gzahn@wm.edu
#SBATCH -o /sciclone/home/gzahn/logs/NEON_SINTAX_%j.out
#SBATCH -e /sciclone/home/gzahn/logs/NEON_SINTAX_%j.err


# Environment setup
NTHREADS=20

# path to vsearch binary
VSEARCH="/sciclone/home/gzahn/programs/vsearch-2.30.1-linux-x86_64/bin/vsearch"
# path to sintax databases
SINTAX_FUNG="/sciclone/scr10/gzahn/DATABASES/Eukaryome/SINTAX_EUK_ITS_v2.0.fasta.gz"
SINTAX_BACT="/sciclone/scr10/gzahn/DATABASES/Silva/silva_138_2_sintax.fa.gz"


# Run vsearch tax assignments
$VSEARCH --sintax ./output/bact_asv.fasta \
         --db $SINTAX_BACT \
         --tabbedout ./output/bact.sintax.tsv \
         --sintax_cutoff 0.8 \
         --strand both \
         --threads $NTHREADS

$VSEARCH --sintax ./output/fung_asv.fasta \
         --db $SINTAX_FUNG \
         --tabbedout ./output/fung.sintax.tsv \
         --sintax_cutoff 0.8 \
         --strand both \
         --threads $NTHREADS

# Run R script to add tax tables from SINTAX to physeq objects
Rscript ./assign_taxonomy_sintax.R


########### Below is backup code for converting newer Silva releases to SINTAX format ##########

# # Converting dada2 silva db to sintax format
# # Input: SILVA DADA2 FASTA
# in="silva_nr99_v138.2_toSpecies_trainset.fa.gz" 
# out="silva_138_2_sintax.fa"
# 
# gzip -cd "$in" | \
# awk '
# BEGIN{
#   FS=";";
#   ranks[1]="d"; ranks[2]="p"; ranks[3]="c";
#   ranks[4]="o"; ranks[5]="f"; ranks[6]="g"; ranks[7]="s";
# }
# # Header line
# /^>/{
#   ++n;
#   # strip initial ">"
#   hdr = substr($0, 2);
#   split(hdr, a, FS);
#   tax = "";
#   for(i=1;i<=7;i++){
#     if(a[i] != "" && a[i] != NULL){
#       # trim whitespace
#       gsub(/^[ \t]+|[ \t]+$/, "", a[i]);
#       if(a[i] != ""){
#         if(tax == "") tax=";tax=";
#         else tax = tax ",";
#         tax = tax ranks[i] ":" a[i];
#       }
#     }
#   }
#   # make a simple unique ID since DADA2 headers lack an accession
#   print ">" "silva1382_ref_" n tax;
#   next;
# }
# # Sequence lines pass through
# { print }
# ' > "$out"




