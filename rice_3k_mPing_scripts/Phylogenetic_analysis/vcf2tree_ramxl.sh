#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=40G
#SBATCH --time=20:00:00
#SBATCH --output=run_qsub.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./


module load samtools
module load vcftools
module load RAxML
vcf=Wildrice_Outgroup_SNP.gatk.snp.pass.rufi48.recode.Hom.norepeat.recode.vcf
output=Wildrice_Outgroup_rufi48_Hom_norepeat_SNP_qsub.raxml

start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=$1
fi

echo "CPU: $CPU"
echo "N: $N"


if [ ! -e $vcf\.tab ]; then
   echo "convert vcf to tab"
   vcf-to-tab < $vcf > $vcf\.tab
fi

if [ ! -e $vcf\.tab.fasta ]; then
   echo "extract fasta from tab"
   perl vcf_tab_to_fasta_alignmentv1.pl --exclude_het -i $vcf\.tab > $vcf\.tab.pos
   perl ~/BigData/software/bin/fastaDeal.pl --attr id:len $vcf\.tab.fasta > $vcf\.tab.fasta.len
fi
   
if [ ! -e Wildrice_Outgroup_Hom_CDS_SNP.raxml ]; then
   python fasta2nexus.py --input $vcf\.tab.fasta
   python ElConcatenero.py -c -if nexus -of phylip -in $vcf\.tab.fasta.nex
   #mv Wildrice_Outgroup_SNP.phy $vcf\.tab.fasta.nex.phy 
   raxmlHPC-PTHREADS-SSE3 -T $CPU -p 12345 -o ./fastq/glaberrima_hap0,./fastq/glumipatula_hap0,./fastq/meridionalis_hap0 -m GTRGAMMA -s $vcf\.tab.fasta.phy -n $output
   raxmlHPC-PTHREADS-SSE3 -T $CPU -p 12345 -x 12345 -# 100 -f a -o ./fastq/glaberrima_hap0,./fastq/glumipatula_hap0,./fastq/meridionalis_hap0 -m GTRGAMMA -s $vcf\.tab.fasta.phy -n $output\_BS
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
