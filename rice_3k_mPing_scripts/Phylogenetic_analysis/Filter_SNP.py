#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import scipy.stats
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python Filter_Het.py --input dbSNP.clean.vcf --vcf Fairchild.gatk.snp.raw.vcf

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

###source=SelectVariants
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Citrus
#scaffold_1	33023	.	T	C	93645.4	PASS	.	GT:AD:DP:GQ:PL	1/1:0,30:30:90:1043,90,0
def Read_VCF(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                info = re.split(r':', unit[9])
                if info[0] == './.':
                    continue
                print info[0], line
                alleles = re.split(r',', info[1])
                alleles.sort(key=int)
                #REF	ALT	GT	#ALLELE0	#ALLELE1	DP
                data['%s:%s' %(unit[0], unit[1])] = [unit[3], unit[4], info[0], alleles[0], alleles[1], info[2]]
    return data

###source=SelectVariants
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Citrus
#scaffold_1	33023	.	T	C	93645.4	PASS	.	GT:AD:DP:GQ:PL	1/1:0,30:30:90:1043,90,0
def Filter_Het(infile):
    data = defaultdict(str)
    HetSNP_vcf = re.sub(r'.vcf', r'.Hom.vcf', infile)
    ofile = open(HetSNP_vcf, 'w')
    print HetSNP_vcf
    r = re.compile(r'0/1')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                if len(unit[4]) == 1 and not r.search(line): #homozygous SNP
                    print >> ofile, line
            else:
                print >> ofile, line
                pass
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--vcf', help='Illumina SNP VCF file')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.vcf) > 0
    except:
        usage()
        sys.exit(2)
  
    #illumina_SNP = Read_VCF(args.vcf)
    Filter_Het(args.vcf)

if __name__ == '__main__':
    main()

