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
python Split_Fastq2PE.py --input rufipogon_W0180_RelocaTE2.te_reads.fq
Input fastq is a mix of read1 and read2 of PE sequence. We split the fatsq file into read1.fq, read2.fq and unpaired.fq using read name.

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

#repeat	Chr1:38668855..38668857	Left_supporting_reads	ERR068809.5052702,ERR068809.7020963,ERR068809.10628656
def locus_reads_list(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                unit[1] = re.sub(r'\.\.', r'_', unit[1])
                unit[1] = re.sub(r':', r'_', unit[1])
                locus = '%s_%s' %(unit[0], unit[1])
                print line, locus
                if len(unit) < 4:
                    continue
                reads = re.split(r',', unit[3])
                for read in reads:
                    data[locus].append(re.split(r':', read)[0])
    prefix = re.sub(r'.list', r'', infile)
    for locus in data.keys():
        ofile = open('%s.%s.list' %(prefix, locus), 'w')
        print >> ofile, '\n'.join(list(set(data[locus])))
        ofile.close()

def get_fastq_seq_by_list(fastqfile, id_list, prefix):
    ofile = open('%s.fq' %(prefix), 'w')
    for record in SeqIO.parse(fastqfile, "fastq"):
        #print 'id:', record.id
        #print 'seq:', record.seq
        unit = re.split(r':', str(record.id))
        record.id = unit[0]
        if id_list.has_key(re.sub(r'read1.', r'', unit[0])) or id_list.has_key(re.sub(r'read2.', r'', unit[0])):
            SeqIO.write(record, ofile, 'fastq')
    ofile.close()

def read_list(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list')
    parser.add_argument('-f', '--fastq')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.list) > 0
    except:
        usage()
        sys.exit(2)

    id_list = read_list(args.list)
    get_fastq_seq_by_list(args.fastq, id_list, re.sub(r'.list', r'', args.list))

if __name__ == '__main__':
    main()

