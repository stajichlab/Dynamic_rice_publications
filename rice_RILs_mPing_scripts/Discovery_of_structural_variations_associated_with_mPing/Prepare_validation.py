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
python Prepare_validation.py --input RIL230_test_filtered

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#2       28724501        28732700        0.00865927
#4       6662401 6663100 0.0185634
#6       1453501 1454200 0.0175716
#7       5266201 5267000 0.0359095
def read_bed(infile, ril, validfile, drawfile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                sv = 'Chr%s:%s-%s' %(unit[0], unit[1], unit[2])
                print >> drawfile, '%s\t%s' %(sv, re.sub(r'GN', r'', ril))
                print >> validfile, 'Chr%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], ril)
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
   
    skip = ['']
    beds = glob.glob('%s/*.6.readdepth_filter_parent.bed' %(args.input))
    validfile = open('%s.validation_table.txt' %(args.input), 'w')
    drawfile  = open('%s.draw.txt' %(args.input), 'w')
    for bed in sorted(beds):
        ril = re.split(r'\.', os.path.split(bed)[1])[0]
        if not ril in skip:
            read_bed(bed, ril, validfile, drawfile)
    validfile.close()
    drawfile.close()

if __name__ == '__main__':
    main()

