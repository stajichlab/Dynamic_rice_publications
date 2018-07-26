#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=fasttree.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

start=`date +%s`

prefix=base_filtered_v0.7.pruneddata_1M_20kb_0.8.reorder.merged.tab
/rhome/cjinfeng/BigData/software/fasttree/FastTree -noml -nome -nt $prefix\.fasta > $prefix\.tree

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

