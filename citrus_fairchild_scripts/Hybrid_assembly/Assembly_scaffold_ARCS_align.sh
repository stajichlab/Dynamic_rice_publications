#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=Step1_ARCS_align.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

export PATH="/rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/arcs/bin:/rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/longranger/longranger-2.1.3:/rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/links_v1.8.5:$PATH"

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


module load samtools/1.3
module load bwa/0.7.12 
#fastq file from longranger basic
fastq=citrus.barcoded.fastq.gz
#scaffold to link
#ref=ecoli_ass.reform.fa
#ref=hsapiens-8reformat.fa
#ref=Fairchild.reform.fa
ref=Fairchild_10xcontig_graph.reform.fa
prefix=citrus
echo "Prepare CHROMIUM interleaved fastq file"
if [ ! -e $prefix\_CHROMIUM_interleaved.fastq.gz ]; then
   echo "......"
   /usr/bin/unpigz -c $fastq | perl -ne 'chomp;$ct++;$ct=1 if($ct>4);if($ct==1){if(/(\@\S+)\sBX\:Z\:(\S{16})/){$flag=1;$head=$1."_".$2;print "$head\n";}else{$flag=0;}}else{print "$_\n" if($flag);}' > $prefix\_CHROMIUM_interleaved.fastq
   pigz $prefix\_CHROMIUM_interleaved.fastq -p $CPU
fi

echo "Align reads by BWA mem"
if [ ! -e $prefix\.CHROMIUM-sorted.bam ]; then
   echo "......"
   bwa index $ref
   bwa mem -t24 $ref -p $prefix\_CHROMIUM_interleaved.fastq.gz | /opt/linux/centos/7.x/x86_64/pkgs/samtools/1.3/bin/samtools view -Sb - | /opt/linux/centos/7.x/x86_64/pkgs/samtools/1.3/bin/samtools sort -n - -o $prefix\.CHROMIUM-sorted.bam
   ls $prefix\.CHROMIUM-sorted.bam > $prefix\.CHROMIUM-sorted.bam.list  
fi

echo "ARCS and LINKS"
C=5      #-c  Minimum number of mapping read pairs/Index required before creating edge in graph. (default: 5)
R=0.05   #-r  Maximum p-value for H/T assignment and link orientation determination. Lower is more stringent (default: 0.05)
E=30000  #-e  End length (bp) of sequences to consider (default: 30000)
L=5      #-l  minimum number of links (k-mer pairs) to compute scaffold (default -l 5, optional)
A=0.9    #-a  maximum link ratio between two best contig pairs (default -a 0.3, optional)
	 #*higher values lead to least accurate scaffolding*
#Fairchild.reform.fa.scaff_s98_c5_l0_d0_e30000_r0.05_original.gv
gv=$ref\.scaff_s98_c5_l0_d0_e30000_r0.05_original.gv
echo $gv
if [ ! -e $gv ]; then
   echo "......"
   arcs -f $ref -a $prefix\.CHROMIUM-sorted.bam.list -s 98 -c $C -l 0 -d 0 -r $R -e $E -v 1 -m 20-10000
   python makeTSVfile.py $gv links_c5r0.05e30000-l5-a0.9.tigpair_checkpoint.tsv $ref
   LINKS -f $ref -s empty.fof -k 20 -b links_c5r0.05e30000-l5-a0.9 -l $L -t 2 -a $A -x 1
   perl rename_scaffold.pl -f links_c5r0.05e30000-l5-a0.9.scaffolds.fa -o links_c5r0.05e30000-l5-a0.9.scaffolds.rename.fa
   perl ~/BigData/software/bin/sumNxx.pl links_c5r0.05e30000-l5-a0.9.scaffolds.rename.fa > links_c5r0.05e30000-l5-a0.9.scaffolds.rename.fa.NXX
   #example
   #arcs -f $ref -a NA24143_genome_phased_namesorted.bam1.sorted.bam.list -s 98 -c $C -l 0 -d 0 -r $R -e $E -v 1 -m 20-10000
   #python makeTSVfile.py hsapiens-8reformat.fa.scaff_s98_c5_l0_d0_e30000_r0.05_original.gv links_c5r0.05e30000-l5-a0.9.tigpair_checkpoint.tsv hsapiens-8reformat.fa
   #LINKS -f $ref -s empty.fof -k 20 -b links_c5r0.05e30000-l5-a0.9 -l $L -t 2 -a $A -x 1
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
