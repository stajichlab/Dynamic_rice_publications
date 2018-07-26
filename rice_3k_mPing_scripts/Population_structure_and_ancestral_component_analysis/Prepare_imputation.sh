#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=Prepare_imputation.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./


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

JAVA=/usr/bin/java
BEAGLE=/rhome/cjinfeng/Rice/Rice_population_sequence/Rice_3000/analysis/Introgression/Beagle/beagle.10Jun18.811.jar
VCF=base_filtered_v0.7.$N\.vcf
OUT=base_filtered_v0.7.$N\.imp

if [ ! -e $OUT\.vcf ]; then

   echo "Processing $OUT"
   $JAVA -Xmx80g -jar $BEAGLE gt=$VCF out=$OUT
   zcat $OUT\.vcf.gz | awk '/^##source/ { printf("##contig=<ID=1,length=43254223>\n##contig=<ID=2,length=35935650>\n##contig=<ID=3,length=36413322>\n##contig=<ID=4,length=35498859>\n##contig=<ID=5,length=29948561>\n##contig=<ID=6,length=31247197>\n##contig=<ID=7,length=29668722>\n##contig=<ID=8,length=28438544>\n##contig=<ID=9,length=22939726>\n##contig=<ID=10,length=23202714>\n##contig=<ID=11,length=29018029>\n##contig=<ID=12,length=27530375>\n");} {print;}' > $OUT\.vcf

fi
end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
