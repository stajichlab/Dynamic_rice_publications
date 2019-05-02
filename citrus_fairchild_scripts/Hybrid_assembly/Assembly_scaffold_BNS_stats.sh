#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=Step1_BNS_stats.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

module load samtools
PATH=$PATH:~/BigData/software/SVcaller/ROOT/bin/
genome=Fairchild_canu1_3.quiver_round1_pilon.fasta

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

echo "Summary BNX bionano raw data"
perl ~/BigData/00.RD/Assembly/bionano/install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/bnx_stats.pl --min_length_kb 150 ~/BigData/00.RD/Assembly/bionano/input/Fairchild/RawMolecules_citrus.bnx
echo "In-silico digest input FASTA into CMAP"
perl ~/BigData/00.RD/Assembly/bionano/install/Solve_03062017Rel/HybridScaffold/03062017/scripts/fa2cmap_multi_color.pl -i ~/BigData/00.RD/Assembly/bionano/input/Pacbio_assembly/genome_ref.fa -e BssSI 1

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
