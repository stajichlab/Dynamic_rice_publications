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
python Ping_Coverage_matrix.py --input bam_5_rufipogon/ --output ping_coverage_rufipogon.matrix

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=10:00:00,mem=10G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = 1
    return data


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

    samtools = '/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.3/bin/samtools'
    bedtools = '~/BigData/software/bedtools2-2.19.0/bin/bedtools'
    bwa      = '/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin/bwa'
    fq_dir   = 'stowaway_check_candidate_fq'
    os.system('mkdir {}'.format(fq_dir))
    bams  = glob.glob('%s/*.bam' %(args.input))
    id_list = readtable('stowaway_check.status.insertion.list')
    #id_list = readtable('stowaway_check.status.insertion.list.fake')

    cmd = []
    for bam in sorted(bams):
        strain    = os.path.split(bam)[1]
        strain    = re.sub(r'.special_ping_locus.bam', '', strain)
        if id_list.has_key(strain):
            cmd.append('{} fastq {} -0 {}/{}.0.fq -1 {}/{}.1.fq -2 {}/{}.2.fq'.format(samtools, bam, fq_dir, strain, fq_dir, strain, fq_dir, strain))
            cmd.append('cat {}/{}.0.fq {}/{}.1.fq {}/{}.2.fq > {}/{}.fq'.format(fq_dir, strain, fq_dir, strain, fq_dir, strain, fq_dir, strain))
            #cmd.append('python Split_Fastq2PE.py --input {}/{}.fq'.format(fq_dir, strain))
            #bwa mem stowaway.flank1k.fa stowaway_check_candidate_fq/B160.fq | samtools view - -b | samtools sort -o test.sort1.bam -
            cmd.append('{} mem stowaway.flank1k.fa {}/{}.fq | {} view - -b | {} sort -o {}/{}.bam -'.format(bwa, fq_dir, strain, samtools, samtools, fq_dir, strain))
            cmd.append('{} index {}/{}.bam'.format(samtools, fq_dir, strain))
            #print '\n'.join(cmd)
    for c in cmd:
        os.system(c)

if __name__ == '__main__':
    main()

