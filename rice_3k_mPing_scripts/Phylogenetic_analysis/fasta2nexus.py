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
python fastPHASE2nexus.py --input test20kb_hapguess_switch.out


    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid



## CX97
#C C C C C C A T T A C G A A C T C T A C G A T A T G C T T A C A A C C T C G G G C G G C T C T A A A C A G C T C C G T C G G G A C G A A
#C C C C C C A T T A C G A A C T C T A C G A T A T G C T T A C A A C C T C G G G C G G C T C T A A A C A G C T C C G T C G G G A C G A A
def convert_nexus(infile):
    data = defaultdict(lambda : str())
    flag = 0
    strain = ''
    hap    = ''
    strain_n  = 0
    loci_n    = 0
    for record in SeqIO.parse(infile, "fasta"):
        strain_n += 1
        loci_n    = len(str(record.seq))
        data[str(record.id)] = str(record.seq)
    ofile = open('%s.nex' %(infile), 'w')
    print >> ofile, '''#NEXUS
Begin data;
Dimensions ntax=%s nchar=%s;
Format datatype=dna missing=? gap=-;
Matrix
''' %(strain_n, loci_n)
    for strain in sorted(data.keys()):
        for i in [0]:
            print >> ofile, '%s_hap%s    %s' %(strain, i, data[strain])
    print >> ofile, ';\nEnd;'

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

    convert_nexus(args.input)

if __name__ == '__main__':
    main()

