#!/usr/bin/python
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
import multiprocessing as mp

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



def split_fastq_seq(fastqfile):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    barcode_sum=defaultdict(lambda : int())
    #read1_s = re.compile(r'read1')
    #read2_s = re.compile(r'read2')
    #fastq_record = defaultdict(lambda )
    #@E00526:73:HJW57ALXX:8:1106:27813:48863_AAACACCAGCGATATA
    for record in SeqIO.parse(fastqfile, "fastq"):
        #print 'id:', record.id
        #print 'seq:', record.seq
        unit = re.split(r'_', str(record.id))
        name = unit[0]
        barcode = unit[1]
        barcode_sum[barcode] += 1
    ofile_barcode_sum = open('%s.barcode.sum' %(fastqfile), 'w')
    for bc in sorted(barcode_sum.keys()):
        print >> ofile_barcode_sum, '%s\t%s' %(bc, barcode_sum[bc])
    ofile_barcode_sum.close()
'''
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
'''

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

def sum_barcode(fastqfile):
    barcode_sum=defaultdict(lambda : int())
    #@E00526:73:HJW57ALXX:8:1106:27813:48863_AAACACCAGCGATATA
    for record in SeqIO.parse(fastqfile, "fastq"):
        #print 'id:', record.id
        #print 'seq:', record.seq
        unit = re.split(r'_', str(record.id))
        name = unit[0]
        barcode = unit[1]
        barcode_sum[barcode] += 1

    barcode_pass = defaultdict(lambda : int())
    ofile_barcode_sum  = open('%s.barcode.sum' %(fastqfile), 'w')
    ofile_barcode_list = open('%s.barcode.list' %(fastqfile), 'w')
    for bc in sorted(barcode_sum.keys()):
        print >> ofile_barcode_list, '%s\t%s' %(bc, barcode_sum[bc])
        if barcode_sum[bc] >= 50:
            barcode_pass[bc] = barcode_sum[bc]
    print >> ofile_barcode_sum, str(len(barcode_pass.keys()))
    ofile_barcode_sum.close()
    ofile_barcode_list.close()
    #os.system("awk '$2 > 50' %s.barcode.list | wc -l > %s.barcode.sum" %(fastqfile, fastqfile))
    #add something here to extract only these molecular with more than 2X50 reads
    newfastq = '%s.50reads.fastq' %(os.path.splitext(fastqfile)[0])
    ofile_fastq = open(newfastq, 'w')
    for record in SeqIO.parse(fastqfile, "fastq"):
        unit = re.split(r'_', str(record.id))
        name = unit[0]
        barcode = unit[1]
        if barcode_pass.has_key(barcode):
            SeqIO.write(record, ofile_fastq, "fastq")
    ofile_fastq.close()
    os.system('rm %s' %(fastqfile))
    return '%s.barcode.sum' %(fastqfile)

def sum_barcode_helper(args):
    return sum_barcode(*args)

##run function with parameters using multiprocess of #cpu
def multiprocess_pool_sum_barcode(parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(sum_barcode_helper, tuple(parameters))
    collect_list = []
    for x in imap_it:
        #print 'status: %s' %(x)
        collect_list.append(x)
    return collect_list

def sum_barcode_by_split_fq(fq_dir, cpu):
    parameters = []
    for fq in glob.glob('%s/*.f*q' %(fq_dir)):
        parameters.append([fq])
        #print fq
    #run sum barcode
    multiprocess_pool_sum_barcode(parameters, cpu) 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-c', '--cpu', default='24')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    fq_prefix = args.input
    if os.path.splitext(args.input)[1] == '.gz':
        fq_prefix = os.path.splitext(args.input)[0]
    print "fq_prefix: %s" %(fq_prefix)
    if not os.path.exists('%s_split' %(fq_prefix)):
        print "split fq into 100,000 reads chunks: %s" %(fq_prefix)
        os.system('perl scripts/fastq_split.pl -s 100_000 -o %s_split %s' %(fq_prefix, args.input))
    sum_barcode_by_split_fq('%s_split' %(fq_prefix), args.cpu)
    os.system('ls %s_split/*.50reads.fastq > %s_split.50reads.list' %(os.path.abspath(fq_prefix), fq_prefix))
 
if __name__ == '__main__':
    main()

