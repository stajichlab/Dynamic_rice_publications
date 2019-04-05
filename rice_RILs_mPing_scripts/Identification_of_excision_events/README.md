## Identification of mPing/Ping/Pong excision

```shell
cp ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Question11_fix_excision_module/*.sh ./
cp ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Question11_fix_excision_module/*.py ./
cp -R ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Question11_fix_excision_module/lib/ ./
```

+ Analyze read coverage at nonreference mPing/Ping insertions using bam files that were aligned to the reference genome and pseudogenome with mPing/Ping inserted in the reference genome.

```shell
python mPing_Boundary_Coverage_MP.py --bam_ref RILs_ALL_bam_correct_merged --bam_pseudo RILs_ALL_unmapped_mping_bam_RILs_514_mPing_Ping --gff_ref RILs.ALL_mPing_Ping.gff --gff_pseudo MSU_r7.Pseudo_mPing_Ping_514.gff --project output_MSU7_nonreference_514_MP3

arguments:
   --bam: directory that has bam files with reads aligned to the reference genome
   --bam_pseudo: directory that has bam files with reads aligned to the pseudogenome
   --gff: gff file of nonreference mPing/Ping insertions in the reference genome
   --gff: gff file of mPing/Ping insertions (become reference insertions) in the pseudogenome
   --project: prefix for the output files
```

