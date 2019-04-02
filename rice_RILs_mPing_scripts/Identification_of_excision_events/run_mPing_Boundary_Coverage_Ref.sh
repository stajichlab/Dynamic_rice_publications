#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=40:00:00
#SBATCH --output=run_mPing_Boundary_Coverage_Ref.sh.%A_%a.stdout
#SBATCH -p intel,batch
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
#python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam_correct_merged --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam_RILs_514_mPing_Ping --gff_ref ../input/RILs.ALL_mPing_Ping.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing_Ping_514.gff --project output_MSU7_nonreference_514 > log 2>&1
python mPing_Boundary_Coverage_Ref.py --bam_ref RILs_ALL_unmapped_mping_bam_Ref_57_mPing_Ping_Pong --bam_pseudo RILs_ALL_bam_correct_merged/ --gff_ref Parent.Pseudo_mPing_Ping_Pong_57.Ref_Shared.gff --gff_pseudo Parent.ALL.mPing_Ping_Pong.Ref_Shared.gff --project output_MSU7_reference_57 > log1 2>&1

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
