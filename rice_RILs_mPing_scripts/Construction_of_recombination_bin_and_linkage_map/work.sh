echo "prepare parental genotype files"
sbatch step01.parent.sh
echo "genotype RILs"
sbatch step01.genotype.sh
echo "construct recombination bin"
sbatch step02.recombination_bin.sh
echo "QTL mapping"
sbatch step03.QTL.sh
