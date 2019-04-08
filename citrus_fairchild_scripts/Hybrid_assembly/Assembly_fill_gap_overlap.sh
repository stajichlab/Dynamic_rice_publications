#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=Graph_Assembly_chr7_merge.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

start=`date +%s`
export PERL5LIB=/rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/perl_lib/lib/site_perl/5.20.2/:/rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10/:$PERL5LIB

chr=chr7
outdir=fairchild_graph_57G_$chr\_bwa_merge
#mkdir fairchild_graph
#time perl /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10//02-PacbioGAPFilter.pl test.blasr.out 10000 95 5000 1000 > fairchild_graph/02-PacbioGAPFilter.txt
#time perl /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10//03-PacbioGAPLinker.pl fairchild_graph/02-PacbioGAPFilter.txt 10000 1000 > fairchild_graph/03-PacbioGAPLinker.txt
#time perl /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10//construct_graph_final.pl fairchild_graph/03-PacbioGAPLinker.txt fairchild_graph/
#time perl /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10//make_ctg_line.pl fairchild_graph/cluster_ori.txt fairchild_graph/cluster_ori_same_chain.txt
#time perl /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10//make_junction_by_pos.pl fairchild_graph/ctg_pairs.txt fairchild_graph/ctg_ctg_ori.txt fairchild_graph/cluster_ori_same_chain.txt fairchild_graph/cluster_ori_same_chain_pos.txt
#time perl /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10//extract_ctg_infor_for_seq.pl fairchild_graph/cluster_ori_same_chain_pos.txt fairchild_graph/cluster_ori_same_chain_pos_for_seq.txt
time perl /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10//extract_seq_by_pos.pl $outdir/cluster_ori_same_chain_pos_for_seq.txt Fairchild_$chr\.fasta_contig.fasta pool_asm.test.10kb.fasta_contig.fasta Final_Connected.fasta /rhome/cjinfeng/BigData/00.RD/Assembly/10xgenomics/tools/ncomms15324-s10/
time cat Final_Connected.fasta non-connected-ctg.fasta > $outdir/Final_Ref.fasta
perl ~/BigData/software/bin/sumNxx.pl $outdir/Final_Ref.fasta > $outdir/Final_Ref.fasta.NXX
mv non-connected-ctg.fasta $outdir/
mv Final_Connected.fasta $outdir/

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
