#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=200G
#SBATCH --time=40:00:00
#SBATCH --output=pilon.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./


#cd $PBS_O_WORKDIR
#qsub -t 1-7 indexbam_qsub.sh

genome=Fairchild_falconv3_20kb_cov2_all_ctg_quiver_round1.fasta
bam=Fairchild.falconv3_20kb_cov2_all_ctg_quiver_round1.bam
prefix=Fairchild_falconv3_20kb_cov2_all_ctg_quiver_round1
outdir=Fairchild_falconv3_20kb_cov2_all_ctg_quiver_round1_pilon
java -Xmx160G -jar /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/pilon/pilon/pilon-1.20.jar --genome $genome --frags $bam --output $prefix --outdir $outdir --threads $SLURM_NTASKS

echo "Done"
