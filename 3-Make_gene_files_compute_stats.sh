#!/bin/bash

#SBATCH -A snic2021-22-847
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -J Makefiles_ComputePiNPiS
#SBATCH --mail-user
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.12
module load vcftools/0.1.16
module load GATK/4.2.0.0
module load SeqKit/0.15.0
module load R_packages/4.1.1

# chromosome:
arg=$1

# start loop: for each gene:
echo "make fasta and vcf for:" 
while read chrm gene start end strand
do
echo $gene
echo $chrm

# Extract gene's fasta
seqkit grep -p $gene /crex/proj/snic2020-16-182/C_rubella/REF/Crubella_183_v1.0.cds.fa > fasta.fasta

# Extract gene's vcf
vcftools --gzvcf /crex/proj/snic2020-16-182/C_rubella/VCF/scaffold_$arg.recode.vcf --chr $chrm --from-bp $start --to-bp $end --recode --out vcf

# Extract gene's gff 
grep $gene /crex/proj/snic2020-16-182/C_rubella/REF/Crubella_183_v1.0.gene.gff3 > $gene.tmp.gff
grep 'CDS' $gene.tmp.gff > gff.gff
rm $gene.tmp.gff

# launch R script 
Rscript ../4-StdFiles_to_PiNPiS.R

# remove used files
rm fasta.fasta vcf.recode.vcf gff.gff full_fasta.fasta vcf.log

# keep potential errors
if grep -i -q 'error' slurm*.out; then
cp slurm*.out > $gene.error.out
cp /dev/null slurm*.out
else
cp /dev/null slurm*.out
fi

done < ../scaffold_$arg.gene.list.txt
