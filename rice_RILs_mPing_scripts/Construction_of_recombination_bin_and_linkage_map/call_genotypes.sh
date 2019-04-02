#correct most at first
#perl make_genotype_makefile.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.redo_UCR.txt -q batch -l "nodes=1:ppn=24,mem=100g,walltime=50:00:00"
#perl make_genotype_makefile_cornell.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.redo_cornell.txt -q batch -l "nodes=1:ppn=24,mem=100g,walltime=50:00:00"

#switch three pairs found by lulu and jinfeng using ping and mping PCR, 89 and 90, 99 and 40, 118 and 119
#perl make_genotype_makefile.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.redo_UCR.switch_20151106.txt -q batch -l "nodes=1:ppn=24,mem=100g,walltime=50:00:00"
#perl make_genotype_makefile_cornell.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.redo_cornell.switch_20151106.txt -q batch -l "nodes=1:ppn=24,mem=100g,walltime=50:00:00"

#new sequence from NextSeq flowcell365, 12 RILs
#genotype_sample_illuminaID.FC365.txt
#perl make_genotype_makefile.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.FC365.txt -q highmem -l "nodes=1:ppn=8,mem=60g,walltime=50:00:00"
#genotype_sample_illuminaID.FC381.txt
#genotype_sample_illuminaID.FC382.txt
#perl make_genotype_makefile.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.FC381.txt -q highmem -l "nodes=1:ppn=8,mem=60g,walltime=50:00:00"
#perl make_genotype_makefile.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.FC382.txt -q highmem -l "nodes=1:ppn=8,mem=60g,walltime=50:00:00"
#perl make_genotype_makefile_gz.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.FC402.txt -q highmem -l "nodes=1:ppn=8,mem=60g,walltime=50:00:00"
##perl make_genotype_makefile_gz.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.FC405.txt -q highmem -l "nodes=1:ppn=8,mem=60g,walltime=50:00:00" #error in genotype_sample_illuminaID.FC405.txt
perl make_genotype_makefile_gz_slurm.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.FC405.txt

#Single Ping RIL from lulu
#perl make_genotype_makefile_gz.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.FC284.txt -q js -l "nodes=1:ppn=8,mem=60g,walltime=50:00:00"
#perl make_genotype_makefile_gz.pl -r /rhome/cjinfeng/Rice/RIL/Illumina_correct -o /rhome/cjinfeng/Rice/RIL/genotypes_correct -f . -i genotype_sample_illuminaID.FC322.txt -q js -l "nodes=1:ppn=8,mem=60g,walltime=50:00:00"
