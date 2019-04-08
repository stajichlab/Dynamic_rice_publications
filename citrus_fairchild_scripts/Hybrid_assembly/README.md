
![Workflow](Hybrid_assembly.pdf)

+ Scaffold contigs with 10x genomics, bionano, and genetic map

```


echo "Genetic map"
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ALLMAPS/Citrus_falconv3_20kb_cov2_p_ctg_quiver_round1_haplomerge_10x_bionano_hybrid
sbatch Assembly_scaffold_genetic_map.sh 

```

+ Fill gaps in the assembly with 10x genomics-PacBio contigs 

```shell
echo "Assemble pooled PacBio long reads isolated by 10x genomics linked reads"
#/rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/bin/assembly_10x_molecular
sbatch Assembly_fill_gap_asm_pooled_reads.sh 

echo "Fill gaps in the assembly (include some manual curation)"
#/rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/bin/contig_gragh/fairchild_graph_57G_chr7_bwa_merge
sbatch Assembly_fill_gap_overlap.sh

echo "Finish gap filling and rename scaffold/contig using chromosomes"
#Working directory: /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/bin/contig_gragh/fairchild_graph_57G_merge_all
sbatch Assembly_fill_gap_finish.sh
```

+ Polish the assembly with pilon for a third round after filling gaps with 10x genomics-PacBio contigs.

```shell
echo "Map illumina paired-end reads of Fairchild to the assembly"
#Working directory: /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/pilon/mapping/
sbatch Assembly_pilon_round3_mapping.sh

echo "Polish the assembly with illumina paired-end reads"
#Working directory: /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/pilon/Citrus/Fairchild_contig_graph.raw.named/
sbatch Assembly_pilon_round3_polish.sh
```
