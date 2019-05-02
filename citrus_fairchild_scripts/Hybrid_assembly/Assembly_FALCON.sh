#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=20G
#SBATCH --time=50:00:00
#SBATCH --output=FALCON_Ecoli.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

start=`date +%s`

module load java/8u25
unset PYTHONPATH
export FALCON_WORKSPACE=/bigdata/stajichlab/cjinfeng/00.RD/Assembly/Pacbio/FALCON_v3.0/virtualenv/FALCON_v3.0/FALCON-integrate
export PYTHONUSERBASE=$FALCON_WORKSPACE/fc_env
export FALCON_PREFIX=${PYTHONUSERBASE}
export PATH=${PYTHONUSERBASE}/bin:${FALCON_PREFIX}/bin:${PATH}

#which fc_run
#fc_run fc_run_ecoli.cfg
fc_run fc_run_large.cfg

#ref=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Ecoli/selfSampleData/reference.fasta
#asm=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/PBcR/K12/9-terminator/asm.ctg.fasta
#/usr/local/bin/dnadiff $ref $asm

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

