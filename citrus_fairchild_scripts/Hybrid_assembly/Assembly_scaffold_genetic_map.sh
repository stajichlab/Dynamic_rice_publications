#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=50G
#SBATCH --time=40:00:00
#SBATCH --output=scaf.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./



start=`date +%s`
module load kent

python -m jcvi.assembly.allmaps merge hk_map.csv np_map.csv lm_map.csv -w weights.txt
python -m jcvi.assembly.allmaps path out.bed scaffolds.fasta -w weights.txt
perl ~/BigData/software/bin/sumNxx.pl out.chr.fasta > out.chr.fasta.NXX
perl ~/BigData/software/bin/sumNxx.pl out.fasta > out.fasta.NXX

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

