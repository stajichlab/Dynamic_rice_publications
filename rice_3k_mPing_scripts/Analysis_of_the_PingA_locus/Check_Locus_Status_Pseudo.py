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
from excision import bamcheck, bamcheck_simple, bamcheck_ref

def usage():
    test="name"
    message='''
python Ping_Coverage_matrix.py --input bam_5_rufipogon/ --output ping_coverage_rufipogon.matrix

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=10:00:00,mem=10G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr6_23521641_23526981_fwd      1       0
#Chr6_23521641_23526981_fwd      2       0
#Chr6_23521641_23526981_fwd      3       0
#Chr6_23521641_23526981_fwd      4       0
def read_cov(infile):
    name     = os.path.split(infile)[1]
    coverage = 0
    data     = []
    count   = 0
    covered = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if int(unit[2]) > 0:
                    data.append('10')
                else:
                    data.append('0')
                count += 1
                if int(unit[2]) > 0:
                    covered += 1
    coverage = covered
    return coverage, data


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

    if not args.output:
        args.output = 'ping_coverage.matrix'

    bams  = glob.glob('%s/*.bam' %(args.input))
    mping = "Chr1_2640500_2640502_fwd:1000-1770"
    bamcheck_file_ref = 'temp_bam_check_file.txt'
    ofile = open('{}.status'.format(args.input), 'w')
    for bam in sorted(bams):
        strain    = os.path.split(bam)[1]
        strain    = re.sub(r'.bam', '', strain)
        ref_flag  = bamcheck_ref(bam, mping, bamcheck_file_ref, strain)
        #1 is insertion, 0 and 4 are excision or no insertion, 2 and 3 are not sure
        print >> ofile, '{}\t{}'.format(strain, ref_flag)
    ofile.close()

if __name__ == '__main__':
    main()

