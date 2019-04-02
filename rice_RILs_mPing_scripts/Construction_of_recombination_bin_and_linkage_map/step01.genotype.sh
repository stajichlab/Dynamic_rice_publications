#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=step01.genotype.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

scripts=`pwd`/scripts

#Trait
perl scripts/trait/subtrait.pl --trait ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL --maqlist BWA.sampleRIL.list > ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.first

#SNP
python $scripts/genotype/Tab2SNP.py --input ../input/fastq/RILs_ALL_bam 
perl $scripts/genotype/RIL_SNP_BWA.pl --fastq ../input/fastq/RILs_ALL_bam --parents NB.RILs.dbSNP.SNPs.parents

echo "done"

