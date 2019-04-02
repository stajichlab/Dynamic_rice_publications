#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=40:00:00
#SBATCH --output=step02.recombination_bin.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

scripts=`pwd`/scripts


export R_LIBS="/bigdata/stajichlab/cjinfeng/software/R_package_MPR/"

#two hours to run for RIL272

#construct recombination bin using MPR package, 3-4 hours
cat $scripts/recombination_map/MPR_hmmrun.R | /opt/linux/centos/7.x/x86_64/pkgs/R/3.2.0/bin/R --slave
#only run R/qtl steps in the scripts to fill and uniq recombination bin, and output to files. Some job was crashed you need to run this.
#deletion the first empty element in MPR.geno.bin!!!!!!!
#cat $scripts/recombination_map/MPR_hmmrun_Rqtl.R | /opt/linux/centos/7.x/x86_64/pkgs/R/3.2.0/bin/R --slave
#draw bin map for each RILs and for each chromosome
#perl $scripts/recombination_map/RIL_drawbin.pl --MPR ./ --chrlen ../input/reference/MSU7.chr.inf

echo "Done"
