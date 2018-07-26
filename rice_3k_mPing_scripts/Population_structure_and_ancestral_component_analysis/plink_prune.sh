#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=plink_prune.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#prefix=3K_coreSNP-v2.1
#prefix=core_v0.7
prefix=base_filtered_v0.7

plink=/rhome/cjinfeng/BigData/software/plink/plink_linux_x86_64_v1.09/plink

#prefix.pruneddata_1M_20kb_0.8
$plink --bfile $prefix --indep-pairwise 1000kb 20kb 0.8
$plink --bfile $prefix --extract plink.prune.in --remove core_v0.7.strain.only.list --make-bed --out $prefix.pruneddata_1M_20kb_0.8

echo "Done"

