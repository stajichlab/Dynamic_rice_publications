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

def usage():
    test="name"
    message='''
python scripts/pool_reads_asm.py --input test.fq_split

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

def write_slurm_shell_canu(fa):
    shell = '%s.canu.sh' %(fa)
    cmd='''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --output=%s.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#run assembly with canu: all in 1
module load java/8u25
canu=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/Canu/canu-1.3/Linux-amd64/bin/canu

$canu -assemble \\
 maxMemory=190 \\
 maxThreads=24 \\
 useGrid=false \\
 -p citrus -d %s_canu \\
 genomeSize=300m \\
 -pacbio-corrected %s

echo "Done"
''' %(shell, fa, fa)

    ofile = open(shell, 'w')
    print >> ofile, cmd
    ofile.close()
    if not os.path.exists('%s_canu' %(fa)):
        job = 'sbatch %s.canu.sh' %(fa)
        print job
        os.system(job)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-c', '--cpu', default='24')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    
    if args.input:
        args.input = os.path.abspath(args.input)   

    ofile_asm_merge = open('%s_asm_merge.fasta' %(args.input), 'w')
    for fa in sorted(glob.glob('%s/*.pacbio_reads.fa_canu/citrus.contigs.fasta' %(args.input))):
        prefix = re.split('\.' , os.path.split(os.path.split(fa)[0])[1])[0]
        print prefix
        for record in SeqIO.parse(fa, "fasta"):
             #print 'id:', record.id
             #print 'seq:', record.seq
             record.id = '%s_%s_%s' %(os.path.split(args.input)[1], prefix, str(record.id))
             SeqIO.write(record, ofile_asm_merge, 'fasta') 
    ofile_asm_merge.close()
    os.system('perl ~/BigData/software/bin/sumNxx.pl %s_asm_merge.fasta > %s_asm_merge.fasta.NXX' %(args.input, args.input))

if __name__ == '__main__':
    main()

