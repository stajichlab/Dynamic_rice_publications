#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=run_speedseq_qsub.sh.%A_%a.stdout
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

module load bcftools

RFMix=/rhome/cjinfeng/Rice/Rice_population_sequence/Rice_3000/analysis/Introgression/RFMix/RFMix_v2/rfmix/rfmix
QRY=base_filtered_v0.7.$N\.imp.QRY.vcf
REF=base_filtered_v0.7.$N\.imp.REF.vcf

bcftools view -S 3k.sample_map.id base_filtered_v0.7.$N\.imp.vcf > base_filtered_v0.7.$N\.imp.REF.vcf
bcftools view -S Query.id base_filtered_v0.7.$N\.imp.vcf > base_filtered_v0.7.$N\.imp.QRY.vcf
$RFMix -f $QRY -r $REF -m 3k.sample_map.txt -g 3k.genetic_map.ALL.txt -o RFMix.3k_$N --chromosome=$N

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
