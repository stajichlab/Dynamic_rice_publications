echo "Genotypeing landrace"
#~/Rice/Rice_population_sequence/Rice_3000/analysis/Genotype_Landrace_SNP270k
bash process_reads_genotype.sh MSU_r7.fa EG4_2 base_filtered_v0.7.pruneddata_1M_20kb_0.8.reorder.vcf ~/Rice/Rice_population_sequence/Rice_3000/analysis/Genotype_Landrace_SNP270k/bam/
bash process_reads_genotype.sh MSU_r7.fa HEG4_2 base_filtered_v0.7.pruneddata_1M_20kb_0.8.reorder.vcf ~/Rice/Rice_population_sequence/Rice_3000/analysis/Genotype_Landrace_SNP270k/bam/
bash process_reads_genotype.sh MSU_r7.fa A123_2 base_filtered_v0.7.pruneddata_1M_20kb_0.8.reorder.vcf ~/Rice/Rice_population_sequence/Rice_3000/analysis/Genotype_Landrace_SNP270k/bam/
bash process_reads_genotype.sh MSU_r7.fa A119_2 base_filtered_v0.7.pruneddata_1M_20kb_0.8.reorder.vcf ~/Rice/Rice_population_sequence/Rice_3000/analysis/Genotype_Landrace_SNP270k/bam/
sbatch vcf_merge.sh
sbatch fasttree.sh


echo "Genotyping O. rufipogon"
sbatch run_speedseq_qsub.sh
python Filter_SNP.py --vcf Wildrice_Outgroup_SNP.gatk.snp.pass.vcf
python Conserved_Gene_SNP.py --gene_list colinear.lst --gtf MSU7.gene.gtf --vcf Wildrice_Outgroup_SNP.gatk.snp.pass.Hom.vcf --repeat_bed /rhome/cjinfeng/Rice/Rice_population_sequence/Rufipogon/reference/MSU_r7.fa.RepeatMasker.out.bed

echo "Phylogenetic analysis of O. rufipogon"
sbatch vcf2tree_ramxl.sh
