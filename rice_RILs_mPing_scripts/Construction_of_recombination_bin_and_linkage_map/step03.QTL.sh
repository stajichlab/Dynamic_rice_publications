#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=step03.QTL.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

scripts=`pwd`/scripts

export R_LIBS="/bigdata/stajichlab/cjinfeng/software/R_package_MPR/"
#export R_LIBS="/rhome/cjinfeng/software/tools/R-2.15.3/library/"

#four hours to run for RIL272

start=`date +%s`

perl $scripts/QTL/RIL_QTL.pl --qtlcart MPR.cross.uniq
#perl $scripts/QTL/RIL_QTL.pl --qtlcart MPR.cross.uniq.test
#cat MPR.cross.uniq.mPing.R | /opt/linux/centos/7.x/x86_64/pkgs/R/3.2.0/bin/R --vanilla --slave


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"


echo "Done"

