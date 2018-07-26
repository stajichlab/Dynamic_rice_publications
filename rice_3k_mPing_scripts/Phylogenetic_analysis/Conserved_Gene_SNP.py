#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python Conserved_Gene_SNP.py --gene_list colinear.lst --gtf MSU7.gene.gtf --vcf Wildrice_Outgroup_SNP.gatk.snp.pass.Hom.vcf
1. Get 1 to 1 colinear gene pairs in BED file;
2. Use this BED to select SNP from VCF

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Os_chr10:
#Os10g01044.1	Ob10g10010.1	1
#Os10g01060.1	Ob10g10030.1	1
#Os10g01080.1	Ob10g10050.1	1
#Os10g01100.1	Ob10g10070.1	1
#Os10g01110.1	Ob10g10080.1	2
def read_gene_list(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                if int(unit[2]) == 1:
                    gene = 'LOC_%s' %(re.split(r'\.', unit[0])[0])
                    data[gene] = 1
    return data

#Chr1	MSU_osa1r7	exon	2903	3268	.	+	.	gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
#Chr1	MSU_osa1r7	exon	3354	3616	.	+	.	gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
def read_gtf(infile, converved_gene):
    data = defaultdict(lambda : int())
    conserved_bed = re.sub(r'.gtf', r'.conserved.bed', infile)
    print conserved_bed
    ofile = open(conserved_bed, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                if unit[2] == 'CDS':
                    gene_id = re.split(r';', unit[8])[0]
                    gene    = re.split(r' ', gene_id)[1]
                    gene    = re.sub(r'"', r'', gene)
                    if converved_gene.has_key(gene):
                        print >> ofile, '%s\t%s\t%s\t%s\t%s' %(unit[0], unit[3], unit[4], unit[2], gene)
    ofile.close()
    return conserved_bed



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene_list')
    parser.add_argument('--gtf')
    parser.add_argument('--vcf')
    parser.add_argument('--repeat_bed')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.vcf) > 0
    except:
        usage()
        sys.exit(2)

    conserved_gene = read_gene_list(args.gene_list)
    conserved_bed  = read_gtf(args.gtf, conserved_gene)
    vcftools = '/opt/linux/centos/7.x/x86_64/pkgs/vcftools/0.1.13/bin/vcftools'
    conserved_vcf = '%s.conserved' %(os.path.splitext(args.vcf)[0])
    norepeat_vcf  = '%s.norepeat' %(os.path.splitext(args.vcf)[0])
    print conserved_vcf
    print norepeat_vcf
    os.system('%s --bed %s --max-missing 0.7 --recode --vcf %s --out %s' %(vcftools, conserved_bed, args.vcf, conserved_vcf))
    os.system('%s --exclude-bed %s --max-missing 0.7 --recode --vcf %s --out %s' %(vcftools, args.repeat_bed, args.vcf, norepeat_vcf))

if __name__ == '__main__':
    main()

