#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=run_speedseq_qsub.sh.stdout
#SBATCH -p batch
#SBATCH --workdir=./


#sbatch --array 1-2 run_speedseq_qsub.sh

module load samtools
genome=/rhome/cjinfeng/Rice/Rice_population_sequence/Rufipogon/reference/MSU_r7.fa

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=$1
fi

FILE=`ls *_1.fastq.gz | grep _1\.fastq\.gz | head -n $N | tail -n 1`
R1=$FILE
R2=`echo $R1 | perl -p -e 's/_1\.fastq/_2.fastq/'`
SAMPLE=${FILE%_1.fastq.gz}

if [ ! -e $SAMPLE\.bam ]; then
echo "mapping $SAMPLE ..."
/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq_20161206/speedseq/bin/speedseq align \
     -t $CPU \
     -o $SAMPLE \
     -R "@RG\tID:id\tSM:$SAMPLE\tLB:lib" \
     $genome \
     $SAMPLE\_1.fastq.gz \
     $SAMPLE\_2.fastq.gz
fi

echo "Done"
