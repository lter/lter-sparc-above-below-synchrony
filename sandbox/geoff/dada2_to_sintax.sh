#!/bin/bash
# Converting dada2 silva db to sintax format
# Input: SILVA DADA2 FASTA
in="silva_nr99_v138.2_toSpecies_trainset.fa.gz"
out="silva_138_2_sintax.fa"

gzip -cd "$in" | \
awk '
BEGIN{
  FS=";";
  ranks[1]="d"; ranks[2]="p"; ranks[3]="c";
  ranks[4]="o"; ranks[5]="f"; ranks[6]="g"; ranks[7]="s";
}
# Header line
/^>/{
  ++n;
  # strip initial ">"
  hdr = substr($0, 2);
  split(hdr, a, FS);
  tax = "";
  for(i=1;i<=7;i++){
    if(a[i] != "" && a[i] != NULL){
      # trim whitespace
      gsub(/^[ \t]+|[ \t]+$/, "", a[i]);
      if(a[i] != ""){
        if(tax == "") tax=";tax=";
        else tax = tax ",";
        tax = tax ranks[i] ":" a[i];
      }
    }
  }
  # make a simple unique ID since DADA2 headers lack an accession
  print ">" "silva1382_ref_" n tax;
  next;
}
# Sequence lines pass through
{ print }
' > "$out"


