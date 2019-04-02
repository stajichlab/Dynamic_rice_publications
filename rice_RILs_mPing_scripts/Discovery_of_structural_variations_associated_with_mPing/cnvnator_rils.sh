#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=50G
#SBATCH --time=100:00:00
#SBATCH --output=stdout
#SBATCH -p intel
#SBATCH --workdir=./

module load samtools
PATH=$PATH:~/BigData/software/SVcaller/ROOT/bin/
start=`date +%s`

for i in `ls ../input/RILs_ALL_bam/*.bam | sed 's/@//'`
do
   echo $i
   prefix=`basename $i .bam`
   prefix1=`echo "${i%.*}"`
   readn=`head -n 1 $prefix1.bam.flagstat | cut -d" " -f1`
   #14880000, 4X
   #18600000, 5X
   if [ ! -e $prefix.readdepth.bed ] && [ $readn -gt 14880000 ]; then
   #echo $prefix
   #echo $prefix1
   #echo $readn
   #samtools view -h $i | sed 's/Chr//g' | samtools view -Sb - -o $prefix.bam
   #samtools index $prefix.bam
   ln -s $i $prefix.bam
   ln -s $i.bai $prefix.bam.bai
   #samtools index $prefix.bam 
   /usr/bin/python2.7 /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator_wrapper.py --cnvnator /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator-multi -T $prefix-cnvnator-temp -t 12 -w 100 -b $prefix.bam -o $prefix.readdepth -c /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//annotations/cnvnator_chroms
   #rm $prefix.bam $prefix.bam.bai
   fi
done


#workdir=/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/lumpy/bin

#/usr/bin/python2.7 /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator_wrapper.py --cnvnator /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator-multi -T $workdir/example/cnvnator-temp -t 12 -w 100 -b $workdir/sample.bam -o $workdir/example/sample.bam.readdepth -c /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//annotations/cnvnator_chroms

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
