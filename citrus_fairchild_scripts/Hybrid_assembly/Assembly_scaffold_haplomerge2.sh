#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=100:00:00
#SBATCH --output=stdout
#SBATCH -p intel
#SBATCH --workdir=./

start=`date +%s`

#### ===========================================================================
#### set the path to executables
#### ===========================================================================
#PATH=/path/to/perl/bin/:/path/to/chainNet/:/path/to/lastz/:/path/to/gapCloser/:/path/to/sspace/:$PATH
PATH=$PATH:/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/HaploMerger2/HaploMerger2_20151124/bin:/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/HaploMerger2/HaploMerger2_20151124/chainNet_jksrc20100603_centOS6:/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/HaploMerger2/HaploMerger2_20151124/lastz_1.02.00_centOS6:/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/SSPACE-LongRead_v1-1/
echo $PATH
export PERL5LIB=""
export PERL5LIB=/rhome/cjinfeng/BigData/software/Perl_lib/lib/site_perl/:$PERL5LIB

#### ===========================================================================
#### lift the file handle limit (require root privilege)
#### ===========================================================================
#ulimit -n 655350

#### ===========================================================================
#### step 1: break the misjoins and output the new diploid assembly
#### ===========================================================================
#### default input:        ./{assembly_name}.fa.gz
#### default final output: ./{assembly_name}_mp.fa.gz
 
./hm.batchA1.initiation_and_all_lastz  genome
./hm.batchA2.chainNet_and_netToMaf     genome
./hm.batchA3.misjoin_processing        genome

#### ===========================================================================
#### step 2: separate/merge two alleles and output haploid assemblies
#### ===========================================================================
#### default input:        ./{assembly_name}_mp.fa.gz
#### default final output: ./{assembly_name}_mp_ref.fa.gz
#### default final output: ./{assembly_name}_mp_alt.fa.gz
  
./hm.batchB1.initiation_and_all_lastz             genome
./hm.batchB2.chainNet_and_netToMaf                genome
./hm.batchB3.haplomerger                          genome
./hm.batchB4.refine_unpaired_sequences            genome
./hm.batchB5.merge_paired_and_unpaired_sequences  genome

#gunzip -d genome_ref.fa.gz
#perl /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/SSPACE-LongRead_v1-1/SSPACE-LongRead.pl -c genome_ref.fa -p preads_norm.fasta -t $PBS_NP -b genome_ref_SSPACELR -a 1000 -i 90 -s 1

#### ===========================================================================
#### step 3: further scaffold the obtained haploid assemblies
#### ===========================================================================
#### default input:        ./{assembly_name}_mp_ref.fa.gz
#### default final output: ./{assembly_name}_mp_ref_ss.fa.gz
#### default final output: ./{assembly_name}_mp_alt_ss.fa.gz
  
#./hm.batchC1.hierarchical_scaffolding                      genome_mp_ref
#./hm.batchC2.update_reference_and_alternative_assemblies   genome_mp_ref

#### ===========================================================================
#### step 4: remove tandem errors from haploid assemblies
#### ===========================================================================
#### default input:        ./{assembly_name}_mp_ref_ss.fa.gz
#### default final output: ./{assembly_name}_mp_ref/alt_ss_rt.fa.gz

#./hm.batchD1.initiation_and_all_lastz     genome_mp_ref
#./hm.batchD2.chainNet_and_netToMaf        genome_mp_ref
#./hm.batchD3.remove_tandem_assemblies     genome_mp_ref

#### ===========================================================================
#### step 5: try to close Ngaps in haploid assemblies
#### ===========================================================================
#### default input:        ./{assembly_name}_mp_ref_ss_rt.fa.gz
#### default final output: ./{assembly_name}_mp_ref/alt_ss_rt_gf.fa.gz


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

#./hm.batchE1.wrapper_for_gapCloser_v1.12  genome_mp_ref/alt_ss_rt
