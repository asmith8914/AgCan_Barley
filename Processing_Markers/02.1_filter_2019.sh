#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=0:10:00
#SBATCH --mail-user=asmith75@uoguelph.ca
#SBATCH --mail-type=ALL

# Objective: Filter prelim_2019 SNPs based on missingness.

module load vcftools/0.1.16

dir="/home/asmith75/scratch/scratch_cedar/AgCan_Barley/Data/filt_2019_then_merge"

vcftools --vcf $dir/prelim_2019_indv_filtered.vcf --max-missing 0.2 --recode --recode-INFO-all --out $dir/prelim_2019_indv_loci_filtered
