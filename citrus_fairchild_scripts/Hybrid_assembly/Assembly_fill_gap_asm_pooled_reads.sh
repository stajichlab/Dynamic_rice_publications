#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=pool_asm.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

module load samtools/1.3
module load bwa/0.7.15

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

#~/BigData/software/seqtk-master/seqtk sample -s seed=11 citrus_CHROMIUM_interleaved.fastq.gz 0.3 > citrus_CHROMIUM_interleaved.0.3.fastq
#python Split_Fastq2PE.py --input citrus_CHROMIUM_interleaved.0.3.fastq
#perl scripts/fastq_split.pl -s 500_000 -o citrus_CHROMIUM_interleaved.fastq_split citrus_CHROMIUM_interleaved.fastq.gz
FILE=p21X.fq.gz
PREFIX=${FILE%.gz}
echo $FILE
echo $PREFIX
python scripts/split_10x_fq.py --input $FILE
python scripts/pool_reads_isolation.py --input $PREFIX\_split --ref citrus.trimmedReads.10kb.fasta
python scripts/pool_reads_asm.py --input $PREFIX\_split
python scripts/pool_reads_asm_merge.py --input $PREFIX\_split

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
