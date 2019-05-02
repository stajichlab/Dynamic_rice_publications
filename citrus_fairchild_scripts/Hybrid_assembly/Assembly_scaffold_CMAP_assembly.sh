#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=20:00:00
#SBATCH --output=Step2_CMAP_assembly.sh.stdout
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

asm=/rhome/cjinfeng/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano_haplotype
ref=/rhome/cjinfeng/BigData/00.RD/Assembly/bionano/input/Pacbio_assembly/genome_ref_BSSSI_0kb_0labels.cmap
perl ~/BigData/00.RD/Assembly/bionano/install/Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/AssembleIrysXeonPhi.pl --assembly_dir $asm --genome 350 -p Fairchild -d
#perl ../install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/cmap_stats.pl --input ~/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano/all_flowcells/bnx_merged_errA_q.cmap
#150kb
#perl ../install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/cmap_stats.pl --input ~/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano/default_t_150/contigs/Fairchild_default_t_150_refineFinal1/FAIRCHILD_DEFAULT_T_150_REFINEFINAL1.cmap > FAIRCHILD_DEFAULT_T_150_REFINEFINAL1.cmap.stat
#perl ../install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/cmap_stats.pl --input ~/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano/relaxed_t_150/contigs/Fairchild_relaxed_t_150_refineFinal1/FAIRCHILD_RELAXED_T_150_REFINEFINAL1.cmap > FAIRCHILD_RELAXED_T_150_REFINEFINAL1.cmap.stat
#perl ../install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/cmap_stats.pl --input ~/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano/strict_t_150/contigs/Fairchild_strict_t_150_refineFinal1/FAIRCHILD_STRICT_T_150_REFINEFINAL1.cmap > FAIRCHILD_STRICT_T_150_REFINEFINAL1.cmap.stat
#180kb
#perl ../install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/cmap_stats.pl --input ~/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano/default_t_180/contigs/Fairchild_default_t_180_refineFinal1/FAIRCHILD_DEFAULT_T_180_REFINEFINAL1.cmap > FAIRCHILD_DEFAULT_T_180_REFINEFINAL1.cmap.stat
#perl ../install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/cmap_stats.pl --input ~/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano/relaxed_t_180/contigs/Fairchild_relaxed_t_180_refineFinal1/FAIRCHILD_RELAXED_T_180_REFINEFINAL1.cmap > FAIRCHILD_RELAXED_T_180_REFINEFINAL1.cmap.stat
#perl ../install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/cmap_stats.pl --input ~/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano/strict_t_180/contigs/Fairchild_strict_t_180_refineFinal1/FAIRCHILD_STRICT_T_180_REFINEFINAL1.cmap > FAIRCHILD_STRICT_T_180_REFINEFINAL1.cmap.stat
#100kb
#perl ../install/Irys-scaffolding/KSU_bioinfo_lab/map_tools/cmap_stats.pl --input ~/BigData/00.RD/Assembly/bionano/bin/Fairchild_bionano/default_t_100/contigs/Fairchild_default_t_100_refineFinal1/FAIRCHILD_DEFAULT_T_100_REFINEFINAL1.cmap > FAIRCHILD_DEFAULT_T_100_REFINEFINAL1.cmap.stat


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
