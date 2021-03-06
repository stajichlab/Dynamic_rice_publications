echo "Call mPing, Ping, Pong insertions in rice"
python ReNameSRA_RelocaTEi_mPing.py --input fastq --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa > log 2>&1
python ReNameSRA_RelocaTEi_mPing.py --input fastq --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/ping.fa > log 2>&1
python ReNameSRA_RelocaTEi_mPing.py --input fastq --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/pong.fa > log 2>&1
python ReNameSRA_sum_Ping.py --input fastq_Ping
python ReNameSRA_sum_Pong.py --input fastq_Pong
python ReNameSRA_sum_mPing.py --input fastq_mPing

python ReNameSRA_RelocaTEi_Rice50.py --input rufipogon --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa > log 2>&1 & 
python ReNameSRA_RelocaTEi_Rice50.py --input rufipogon --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/ping.fa > log 2>&1 &
python ReNameSRA_RelocaTEi_Rice50.py --input rufipogon --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/pong.fa > log 2>&1 &
python ReNameSRA_sum_Ping.py --input rufipogon_RelocaTE2_Ping > log 2>&1 &
python ReNameSRA_sum_Pong.py --input rufipogon_RelocaTE2_Pong > log 2>&1 &
python ReNameSRA_sum_mPing.py --input rufipogon_RelocaTE2_mPing > log 2>&1 &

echo "Estimate mPing, Ping, Pong copy numbers in rice and O. rufipogon using a read-depth method"
python Rice3k_copy_number_depth_window_mPing.py --input Rice3k_3000_RelocaTEi_mPing_NM2 --depth ping_coverage_3k.bam.summary
python Rice3k_copy_number_depth_window_Ping.py --input Rice3k_3000_RelocaTEi_Ping_NM2 --depth ping_coverage_3k.bam.summary
python Rice3k_copy_number_depth_window_Ping.py --input Rice3k_3000_RelocaTEi_Pong_NM2 --depth ping_coverage_3k.bam.summary

python Rice3k_copy_number_depth_window_mPing.py --input Rufipogon_57_RelocaTE2_mPing_NM2/ --depth Rufipogon_57.bam.summary
python Rice3k_copy_number_depth_window_Ping.py --input Rufipogon_57_RelocaTE2_Pong_NM2/ --depth Rufipogon_57.bam.summary
python Rice3k_copy_number_depth_window_Ping.py --input Rufipogon_57_RelocaTE2_Ping_NM2/ --depth Rufipogon_57.bam.summary
