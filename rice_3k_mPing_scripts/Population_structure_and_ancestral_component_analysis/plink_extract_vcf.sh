#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=plink_makefasta.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

prefix=base_filtered_v0.7

plink=/rhome/cjinfeng/BigData/software/plink/plink_linux_x86_64_v1.09/plink

if [ ! -e $prefix.vcf ]; then 

echo "extract vcf for individual chromosome"
for chr in {1..12};
do
$plink --bfile $prefix --recode vcf-iid --chr $chr --out $prefix.$chr
done

fi

echo "Done"

