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
python ../Change_Fastq_Header.py --fastq ERS469852_1.te_repeat.ContainingReads.fq --header read1

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

def header_change(fastqfile, prefix):
    ofile = open(re.sub(r'.fq', r'.name.fq', fastqfile), 'w')
    for record in SeqIO.parse(fastqfile, "fastq"):
        #print 'id:', record.id
        #print 'seq:', record.seq
        #unit = re.split(r':', str(record.id))
        record.id = re.split(r':', record.id)[0]
        record.id = '%s.%s' %(prefix, record.id)
        record.description = ''
        SeqIO.write(record, ofile, 'fastq')
    ofile.close()

def header_change2(fastqfile, prefix):
    ofile = open(re.sub(r'.fq', r'.name.fq', fastqfile), 'w')
    line_n= 0
    with open (fastqfile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r' ', line)
            if len(line) > 2:
                line_n += 1
                if line_n == 1:
                    header = re.sub(r'@', r'', unit[0])
                    print >> ofile, '@%s.%s' %(prefix, header)
                elif line_n == 2:
                    print >> ofile, line
                elif line_n == 3:
                    print >> ofile, '+'
                elif line_n == 4:
                    print >> ofile, line
                    line_n = 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--header')
    parser.add_argument('-f', '--fastq')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.fastq) > 0
    except:
        usage()
        sys.exit(2)

    header_change2(args.fastq, args.header)

if __name__ == '__main__':
    main()

