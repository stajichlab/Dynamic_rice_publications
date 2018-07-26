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
python Ping_locus_reads_list.py --input ERS467761_mPing_assembly/ERS467761.repeat.reads.list

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

#Chr3	ERS469922	RelocaTE_i	33213691	33213693 .	-	.	ID=repeat_Chr3_33213691_33213693;Strain=ERS469922;
def readgff(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                index= '%s.repeat_%s_%s_%s' %(unit[1], unit[0], unit[3], unit[4])
                temp = re.split(r';', unit[8])
                index2 = re.sub(r'ID=', r'%s.' %(unit[1]), temp[0])
                if not index == index2:
                    index = index2
                data[index] = 1
    return data

#repeat	Chr1:38668855..38668857	Left_supporting_reads	ERR068809.5052702,ERR068809.7020963,ERR068809.10628656
def locus_reads_list(infile, ping_loci):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                unit[1] = re.sub(r'\.\.', r'_', unit[1])
                unit[1] = re.sub(r':', r'_', unit[1])
                locus = '%s_%s' %(unit[0], unit[1])
                #print line, locus
                if len(unit) < 4:
                    continue
                reads = re.split(r',', unit[3])
                for read in reads:
                    data[locus].append(re.split(r':', read)[0])
    prefix = re.sub(r'.repeat.reads.list', r'', infile)
    for locus in data.keys():
        temp = '%s.%s' %(os.path.split(prefix)[1], locus)
        if not ping_loci.has_key(temp):
            continue
        print temp
        ofile = open('%s.%s.list' %(prefix, locus), 'w')
        print >> ofile, '\n'.join(list(set(data[locus])))
        ofile.close()

def split_fastq_seq(fastqfile):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    read1_s = re.compile(r'read1')
    read2_s = re.compile(r'read2')
    #fastq_record = defaultdict(lambda )
    for record in SeqIO.parse(fastqfile, "fastq"):
        #print 'id:', record.id
        #print 'seq:', record.seq
        unit = re.split(r':', str(record.id))
        if read1_s.match(unit[0]):
            readname = re.sub(r'read1.', r'', unit[0])
            record.id = readname
            record.description = ''
            data[readname][1] = record
            #print 'read1:', readname
        elif read2_s.match(unit[0]):
            readname = re.sub(r'read2.', r'', unit[0])
            record.id = readname
            record.description = ''
            data[readname][2] = record
            #print 'read2:', readname

    read1_fq = re.sub(r'.fq', r'_1.fq', fastqfile)
    read2_fq = re.sub(r'.fq', r'_2.fq', fastqfile)
    unpaired_fq = re.sub(r'.fq', r'_unpaired.fq', fastqfile)
    ofile_1 = open(read1_fq, 'w')
    ofile_2 = open(read2_fq, 'w')
    ofile_unpaired = open(unpaired_fq, 'w')
    for read in sorted(data.keys()):
        if data[read].has_key(1) and data[read].has_key(2):
            #print >> ofile_1, 'paired: %s' %(read)
            #print >> ofile_2, 'paired: %s' %(read)
            SeqIO.write(data[read][1], ofile_1, 'fastq')
            SeqIO.write(data[read][2], ofile_2, 'fastq')
        elif data[read].has_key(1):
            #print >> ofile_unpaired, 'read1: %s' %(read)
            SeqIO.write(data[read][1], ofile_unpaired, 'fastq')
        elif data[read].has_key(2):
            #print >> ofile_unpaired, 'read2: %s' %(read)
            SeqIO.write(data[read][2], ofile_unpaired, 'fastq')
        else:
            print 'error: %s' %(read)
    ofile_1.close()
    ofile_2.close()
    ofile_unpaired.close() 

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
    parser.add_argument('-g', '--gff')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.gff:
        args.gff = 'Rice3k_3000_RelocaTEi_Ping.CombinedGFF.ALL.high_conf.gff'

    ping_loci = readgff(args.gff)
    locus_reads_list(args.input, ping_loci)

if __name__ == '__main__':
    main()

