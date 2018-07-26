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


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


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


    #fq_RelocaTE2_Ping_NM2_PE_assembly/rufipogon_W0106.assembly/contigs.fa
    #NODE_2_length_548_cov_12.135036
    ofile = open('%s.contigs.fa' %(args.input), 'w')
    s = re.compile(r'length_(\d+)_cov')
    contigs = glob.glob('%s/*/*/contigs.fa' %(args.input))
    for ctg in sorted(contigs):
        strain = os.path.split(os.path.split(ctg)[0])[1]
        strain = re.sub(r'.assembly', r'', strain)
        strain = re.sub(r'.repeat.reads.', '_', strain)
        for record in SeqIO.parse(ctg,"fasta"):
            length = s.search(record.id).groups(0)[0] if s.search(record.id) else 0
            record.id = '%s_len%s' %(strain, str(length))
            record.anno = ''
            SeqIO.write(record, ofile, 'fasta')
    ofile.close()

if __name__ == '__main__':
    main()

