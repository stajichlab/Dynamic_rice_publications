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
python CircosConf.py --input circos.config --output pipe.conf

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

##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ERS467817.repeat_Chr9_20454261_20454263 node11
#chrUn   14      .       G       A       .       .       DP=31   GT:DP   0/0:1   0/0:1
def vcf2txt(infile):
    data = defaultdict(lambda : list())
    node = []
    position = []
    genotype = defaultdict(lambda : list())
    node_gt  = defaultdict(lambda : defaultdict(lambda : str()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            #print line
            if len(line) > 2 and line.startswith(r'#CHROM'): 
                unit = re.split(r'\t',line)
                for i in range(9, len(unit)):
                    node.append(unit[i])
            elif len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t',line)
                position.append(unit[1])
                genotype[unit[1]].append(unit[3]) 
                genotype[unit[1]].extend(re.split(r',', unit[4]))
                #print genotype[unit[1]]
                for i in range(9, len(unit)):
                    #node_gt[node[i-9]] = unit[i][0]
                    #print node[i-9], unit[i][0]
                    if unit[i][0] == 0:
                        node_gt[node[i-9]][unit[1]] = genotype[unit[1]][int(unit[i][0])]
                    else:
                        node_gt[node[i-9]][unit[1]] = genotype[unit[1]][int(unit[i][0])]
    ofile = open('%s.haplotype.table' %(infile), 'w')
    print >> ofile, 'Position\t%s' %(','.join(position)) 
    for node_family in node:
        base = []
        for pos in position:
            base.append(node_gt[node_family][pos])
        print >> ofile, '%s\t%s' %(node_family, ','.join(base))
    ofile.close()
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    vcf2txt(args.input)

if __name__ == '__main__':
    main()

