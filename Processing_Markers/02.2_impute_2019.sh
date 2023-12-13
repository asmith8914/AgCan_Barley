#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=0:30:00
#SBATCH --mail-user=asmith75@uoguelph.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=4G

# Objective: Impute filtered prelim_2019 using BEAGLE.

module load StdEnv/2020
module load beagle/5.4

dir="/home/asmith75/scratch/scratch_cedar/AgCan_Barley/Data/filt_2019_then_merge"

java -jar ${EBROOTBEAGLE}/beagle.22Jul22.46e.jar gt=$dir/prelim_2019_indv_loci_filtered.recode.vcf out=$dir/prelim_2019_filt_impute
