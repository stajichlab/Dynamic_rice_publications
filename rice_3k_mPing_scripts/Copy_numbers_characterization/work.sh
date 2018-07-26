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

echo "Estimate mPing, Ping, Pong copy numbers in rice and O. rufipogon"
python Rice3k_copy_number_depth.py --input Rice3k_3000_RelocaTEi_Ping_NM2 --output Rice3k_3000_RelocaTEi_Ping_NM2.Ping_copy.txt &
python merge_table_ping_and_depth_3k_mapped_depth.py --input Rice3k_3000_RelocaTEi_Ping_NM2.Ping_copy.txt --depth rice_3k_depth.txt
python Rice3k_copy_number_depth.py --input Rice3k_3000_RelocaTEi_Pong_NM2 --output Rice3k_3000_RelocaTEi_Pong_NM2.Pong_copy.txt &
python merge_table_ping_and_depth_3k_mapped_depth.py --input Rice3k_3000_RelocaTEi_Pong_NM2.Pong_copy.txt --depth rice_3k_depth.txt
python Rice3k_copy_number_depth_mPing.py --input Rice3k_3000_RelocaTEi_mPing_NM2 --output Rice3k_3000_RelocaTEi_mPing_NM2.mPing_copy.txt &
python merge_table_ping_and_depth_3k_mapped_depth.py --input Rice3k_3000_RelocaTEi_mPing_NM2.mPing_copy.txt --depth rice_3k_depth.txt

python Rice3k_copy_number_depth.py --input Rufipogon_RelocaTE2_Ping_NM2 --output Rufipogon_RelocaTE2_Ping_NM2.Ping_copy.txt &
python merge_table_ping_and_depth.py --input Rufipogon_RelocaTE2_Ping_NM2.Ping_copy.txt --depth Rufipogon.bam.summary
python Rice3k_copy_number_depth.py --input Rufipogon_RelocaTE2_Pong_NM2 --output Rufipogon_RelocaTE2_Pong_NM2.Pong_copy.txt &
python merge_table_ping_and_depth.py --input Rufipogon_RelocaTE2_Pong_NM2.Pong_copy.txt --depth Rufipogon.bam.summary
python Rice3k_copy_number_depth_mPing.py --input Rufipogon_RelocaTE2_mPing_NM2 --output Rufipogon_57_RelocaTE2_mPing_NM2.mPing_copy.txt &
python merge_table_ping_and_depth.py --input Rufipogon_RelocaTE2_mPing_NM2.mPing_copy.txt --depth Rufipogon.bam.summary


