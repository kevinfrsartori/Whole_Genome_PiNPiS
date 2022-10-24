#!/bin/bash

#SBATCH -A snic2021-22-847
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -J Make_filter_vcf
#SBATCH --mail-user 
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load vcftools/0.1.16

# your own filtering options
vcftools --vcf /crex/proj/snic2020-16-182/nobackup/private/marion/full.scaff.geno.rawSNPs.vcf --out /crex/proj/snic2020-16-182/C_rubella/VCF/Carub.full --recode --keep /crex/proj/snic2020-16-182/nobackup/private/marion/lists/listCrComplete.txt
