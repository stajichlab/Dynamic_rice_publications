#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=run_speedseq_qsub.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


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

sed 's/>super-scaffold1/>chr1/' ../fairchild_graph_57G_chr1_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>super-scaffold1/>chr2/' ../fairchild_graph_57G_chr2_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>super-scaffold1/>chr3/' ../fairchild_graph_57G_chr3_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>super-scaffold1/>chr4/' ../fairchild_graph_57G_chr4_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>super-scaffold1/>chr5/' ../fairchild_graph_57G_chr5_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>super-scaffold1/>chr6/' ../fairchild_graph_57G_chr6_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>super-scaffold1/>chr7/' ../fairchild_graph_57G_chr7_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>super-scaffold1/>chr8/' ../fairchild_graph_57G_chr8_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>super-scaffold1/>chr9/' ../fairchild_graph_57G_chr9_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.fasta >> Fairchild_contig_graph.raw.fasta

cp Fairchild_contig_graph.raw.fasta Fairchild_contig_graph.raw_chr.fasta
perl ~/BigData/software/bin/sumNxx.pl Fairchild_contig_graph.raw_chr.fasta > Fairchild_contig_graph.raw_chr.fasta.NXX
sed 's/>/>chrUN_/' ../Fairchild_chrUN.fasta >> Fairchild_contig_graph.raw.fasta
sed 's/>/>chrUN_/' ../fairchild_graph_57G_chr*_bwa_merge/Final_Ref.Self.Update_AGP.table.agp.leftover.fasta >> Fairchild_contig_graph.raw.fasta
perl ~/BigData/software/bin/sumNxx.pl Fairchild_contig_graph.raw.fasta > Fairchild_contig_graph.raw.fasta.NXX
python ChangeName.py --input Fairchild_contig_graph.raw.fasta --output Fairchild_contig_graph.raw.named.fasta
perl ~/BigData/software/bin/fastaDeal.pl --attr id:len Fairchild_contig_graph.raw.named.fasta > Fairchild_contig_graph.raw.named.fasta.len
perl ~/BigData/software/bin/sumNxx.pl Fairchild_contig_graph.raw.named.fasta > Fairchild_contig_graph.raw.named.fasta.NXX

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
