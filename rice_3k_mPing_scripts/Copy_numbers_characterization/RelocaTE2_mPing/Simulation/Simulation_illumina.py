#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Simulation_illumina.py --input ./MSU7.Chr4.mPing --size 500

    '''
    print message

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
    parser.add_argument('-s', '--size')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.size:
        args.size = 500
 
    read_len = 100
    depth  = [1, 2]
    #depth = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]
    pirs = '/rhome/cjinfeng/software/tools/TE_diversity/pirs/pirs'
    files = glob.glob('%s/*.fasta' %(args.input))
    for fa in files:
        print fa
        read_dir = '%s_reads' %(os.path.splitext(fa)[0])
        if not os.path.exists(read_dir):
            os.mkdir(read_dir)
        for dep in depth:
            read_pre = '%s_%s' %(read_dir, dep)
            cmd = '%s simulate %s -m %s -x %s -l %s -o %s > %s_%s_%s.o 2> %s_%s_%s.e' %(pirs, fa, args.size, dep, read_len, read_pre, read_pre, read_len, args.size, read_pre, read_len, args.size)
            mv  = 'mv %s.* %s' %(read_dir, read_dir)
            os.system(cmd)
            os.system(mv)
            #print '%s\n%s' %(cmd, mv)

if __name__ == '__main__':
    main()

