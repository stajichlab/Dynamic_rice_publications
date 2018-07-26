echo "population structure analysis"
sbatch plink_prune.sh
sbatch --array 2-10 run_ADMIXTURE.sh
python Sort_ClusterQ8.py --strain_list core_v0.7.pruneddata3.structure.nosex --meanQ_list base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q | sort -k2,2nr -k4,4n -k5,5n > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.txt
python Summary_Cluster_Color.py --strain_list core_v0.7.pruneddata3.structure.nosex --meanQ_list  base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q
grep "Indica" base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.txt| sort -k4,4n -k5,5n > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.IND.txt
grep "Aus" base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.txt| sort -k4,4nr -k5,5nr > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.AUS.txt
grep "Basmati" base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.txt| sort -k4,4nr -k5,5nr > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.ARO.txt
grep "Temp" base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.txt| sort -k4,4n -k5,5n > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.TEJ.txt
grep "Trop" base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.txt| sort -k4,4nr -k5,5nr > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.TRJ.txt
cat base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.TRJ.txt base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.TEJ.txt base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.ARO.txt base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.AUS.txt base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.IND.txt > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.merge.txt
python Reorder_Q.py --strain_list core_v0.7.pruneddata3.structure.nosex --meanQ_list base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q --reorder_list base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.merge.txt
paste base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.merge.txt base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.reordered.txt | cut -f1,6- > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.reordered.strains.txt
cat 3K_coreSNP-v2.1.binary.tab.landrace.nj.landrace.tree_mPing_Ping_Pong.ADMIXTRUE.R| R --slave
python Assign_Individual_Pop.py --Qtable  base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.reordered.strains.txt --colortable base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.color.txt > base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.subgroup.assign.txt

echo "ancestral analysis"
sbatch plink_extract_vcf.sh
sbatch --array 1-12 Prepare_imputation.sh
sbatch --array 1-12 Prepare_Run.sh
