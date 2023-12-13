#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=0:0:15
#SBATCH --mail-user=asmith75@uoguelph.ca
#SBATCH --mail-type=ALL

# Objective: Convert vcf to 012 format.

module load vcftools/0.1.16

dir="/home/asmith75/scratch/AgCan_Barley"

vcftools --vcf $dir/prelim_merged_commonsites.vcf --012 --out $dir/prelim_merged_comsites_012
