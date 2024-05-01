# Overview

ChIPseq analysis of the human transcription factor Runx1

## Commands for initializing environment and executing the pipeline

conda env create -f base_env.yml
snakemake -s your_snakefile.snake --sdm conda -c 1

snakemake -s your_snakefile.snake --sdm conda --executor cluster-generic \
--cluster-generic-submit-cmd "qsub -P bf528 -pe omp {threads}" --jobs X

## Pipeline files 

The week 1-4 files provide the code for the entire pipeline for the ChiP-seq analysis 
