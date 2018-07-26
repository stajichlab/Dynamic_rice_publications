#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=50G
#SBATCH --time=20:00:00
#SBATCH --output=run_ADMIXTURE.sh.%A_%a.stdout
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

admixture=/rhome/cjinfeng/BigData/software/ADMIXTURE/admixture_linux-1.3.0/admixture
seed=2
$admixture base_filtered_v0.7.pruneddata_1M_20kb_0.8.bed $N -j$CPU -s $seed --cv

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
