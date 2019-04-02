echo "map reads to reference genome"
perl make_genotype_makefile_gz_slurm.pl -r /rhome/cjinfeng/Rice/RIL/Illumina -o /rhome/cjinfeng/Rice/RIL/genotypes -f . -i genotype_sample_illuminaID.txt
echo "prepare parental genotype files"
sbatch step01.parent.sh
echo "genotype RILs"
sbatch step01.genotype.sh
echo "construct recombination bin"
sbatch step02.recombination_bin.sh
echo "QTL mapping"
sbatch step03.QTL.sh
