#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=run_speedseq_falconv3_20kb_cov2_all_ctg_quiver_round1.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#cd $PBS_O_WORKDIR
module load samtools
PATH=$PATH:~/BigData/software/SVcaller/ROOT/bin/

start=`date +%s`
RIL=Fairchild
#genome=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa
genome=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/pilon/mapping/Fairchild_falconv3_20kb_cov2_all_ctg_quiver_round1.fasta

# Example speedseq commands on a small slice of chromosome 20

# 1. Align with BWA
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq align -t 16 -o $RIL -R "@RG\tID:id\tSM:$RIL\tLB:lib"\
#   MSU_r7.fa \
#   $RIL\.1.fq \
#   $RIL\.2.fq

#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq align \
/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq_20161206/speedseq/bin/speedseq align \
     -t $SLURM_NTASKS \
     -o $RIL.falconv3_20kb_cov2_all_ctg_quiver_round1 \
     -R "@RG\tID:id\tSM:$RIL\tLB:lib" \
     $genome \
     $RIL\_1.fastq.gz \
     $RIL\_2.fastq.gz

# 2. Detect SNVs and indels
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq var \
#    -o $RIL \
#    $genome \
#    $RIL.bam

# 3. Detect SVs
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq sv \
#    -o $RIL \
#    -B $RIL\.bam \
#    -S $RIL\.splitters.bam \
#    -D $RIL\.discordants.bam \
#    -R /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/lumpy/bin/MSU_r7.fa \
#    -t $PBS_NP \
#    -d \
#    -k \
#    -v

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
