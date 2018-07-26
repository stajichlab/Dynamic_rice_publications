echo "Locus specific analysis of Ping +16G/A SNP in rice"
python Ping_SNP_PE_locus_reads.py --input Rice3k_3000_RelocaTEi_Ping > log 2>&1

echo "Analysis of Ping +16G/A SNP in O. rufipogon"
python Ping_SNP_PE.py --input Rufipogon_57_RelocaTE2_Ping > log 2>&1
