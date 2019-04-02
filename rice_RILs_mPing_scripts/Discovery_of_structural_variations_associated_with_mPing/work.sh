echo "call SV using cnvnator"
sbatch cnvnator_rils.sh
echo "filter by removing SVs overlaping with black regions on reference genome"
python Filter_CNVnator.py --input RIL272
echo "prepare candidate for manual inspection and validation"
python Prepare_validation.py --input RIL272_filtered
