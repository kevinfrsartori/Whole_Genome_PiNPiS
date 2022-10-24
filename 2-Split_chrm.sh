#!/bin/bash

#SBATCH -A snic2021-22-847
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -J Split_chrm
#SBATCH --mail-user 
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.12
module load vcftools/0.1.16
module load GATK/4.2.0.0


chrm=$1

#make vcf per chromosome, filter out non variant and indels
vcftools --vcf /crex/proj/snic2020-16-182/C_rubella/VCF/Carub.full.recode.vcf --out /crex/proj/snic2020-16-182/C_rubella/VCF/scaffold_$chrm --chr scaffold_$chrm --recode --remove-indels --mac 1
#make index
gatk --java-options "-Xmx8G" IndexFeatureFile -I /crex/proj/snic2020-16-182/C_rubella/VCF/scaffold_$chrm.recode.vcf

