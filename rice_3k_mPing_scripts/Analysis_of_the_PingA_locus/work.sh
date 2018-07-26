echo "Analyze Stowaway element insertion"
python Check_subbam.py --input rice_line_ALL_3000.anno.list
python Check_Locus_Status.py --input stowaway_check > stowaway_check.status
python Check_bam_to_fq.py --input stowaway_check
python Check_Locus_Status_Pseudo.py --input stowaway_check_candidate_fq > stowaway_check_candidate_fq.log

echo "Analyze PingA insertion"
python mPing_Boundary_Coverage_PingA.py --bam_ref 3k_stowaway_strains_Ref_bam --bam_pseudo 3k_stowaway_strains_Pseudo_bam > 3k_stowaway_strains.status

