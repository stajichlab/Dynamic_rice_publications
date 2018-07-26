echo "assemble mPing"
python Ping_SNP_PE_asm.py --input Rice3k_3000_RelocaTEi_Ping > log 2>&1
echo "merge assembled sequences"
python merge_assembly_locus.py --input Rice3k_3000_RelocaTEi_Ping_NM2_PE_Assembly
echo "extract mPing candidate from assembled sequences using TARGeT"
python make_target_nonauto_general_redo.py ~/BigData/00.RD/Transposon_Oryza/OGE_genomes/mPing_Target/query ~/BigData/00.RD/Transposon_Oryza/OGE_genomes/mPing_Target/Target ~/BigData/00.RD/Transposon_Oryza/OGE_genomes/mPing_Target/reference/MSU7.fa Target_Run_mPing_MSU7
echo "Clean and classify assembled mPing sequences"
python Extract_mPing_varients.py --input rice3k_mPing_target

