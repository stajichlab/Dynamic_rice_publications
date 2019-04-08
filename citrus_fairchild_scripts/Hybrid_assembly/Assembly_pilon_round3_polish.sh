#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=200G
#SBATCH --time=4:00:00
#SBATCH --output=pilon.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./


#cd $PBS_O_WORKDIR
#qsub -t 1-7 indexbam_qsub.sh

prefix=Fairchild_contig_graph.raw.named
outdir=$prefix\_pilon
genome=$prefix\.fasta
bam=$prefix\.bam
java -Xmx160G -jar /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/pilon/pilon/pilon-1.20.jar --genome $genome --frags $bam --output $prefix --outdir $outdir --threads $SLURM_NTASKS

echo "Done"
