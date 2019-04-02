#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=10:00:00
#SBATCH --output=step01.parent.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

scripts=`pwd`/scripts

#SNP called by Sofia
perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.VQSR.vcf


echo "done"

