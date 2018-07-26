#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=10G
#SBATCH --time=40:00:00
#SBATCH --output=index_vcf.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./



start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

module load vcftools
vcf-merge base_filtered_v0.7.pruneddata_1M_20kb_0.8.reorder.vcf.gz HEG4_2.genotype.vcf.gz EG4_2.genotype.vcf.gz A123_1.genotype.vcf.gz A119_2.genotype.vcf.gz > base_filtered_v0.7.pruneddata_1M_20kb_0.8.reorder.merged.vcf 

prefix=base_filtered_v0.7.pruneddata_1M_20kb_0.8.reorder.merged
vcf-to-tab < $prefix.vcf > $prefix.tab
perl vcf_tab_to_fasta_alignmentv1.pl -i $prefix.tab > $prefix.tab.pos


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
