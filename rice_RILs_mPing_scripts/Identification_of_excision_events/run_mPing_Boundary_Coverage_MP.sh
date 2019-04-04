#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=run_mPing_Boundary_Coverage_MP.sh.%A_%a.stdout
#SBATCH -p intel,batch,stajichlab
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

#python mPing_Boundary_Coverage.py --bam_ref bam_MSU7 --bam_pseudo bam_MSU7_Pseudo_insertion --gff_ref HEG4.ALL_mPing_Ping.gff --gff_pseudo MSU_r7.Pseudo_mPing_Ping_422.gff --project output_MSU7_nonreference_422 > log 2>&1
python mPing_Boundary_Coverage_MP.py --bam_ref RILs_ALL_bam_correct_merged --bam_pseudo RILs_ALL_unmapped_mping_bam_RILs_514_mPing_Ping --gff_ref RILs.ALL_mPing_Ping.gff --gff_pseudo MSU_r7.Pseudo_mPing_Ping_514.gff --project output_MSU7_nonreference_514_MP3 --cpu $CPU
python Ping_number_RILs.High_exicison.py --csv output_MSU7_nonreference_514_MP3_mPing --ping_code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt --output output_MSU7_nonreference_514_MP3_mPing_GT_Ping_code

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
