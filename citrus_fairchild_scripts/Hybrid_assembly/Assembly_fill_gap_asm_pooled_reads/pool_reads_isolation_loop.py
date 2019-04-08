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
python scripts/pool_reads_isolation.py --input test.fq.gz --ref citrus.trimmedReads.10kb.fasta --cpu 24

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

def write_slurm_shell_bwa(fq, ref, cpu):
    shell = '%s.sh' %(fq)
    cmd='''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=%s
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=%s.stdout
#SBATCH -p intel
#SBATCH --workdir=./

module load bamtools/2.4.0
module load samtools/1.3
module load bwa/0.7.15
module load seqtk

CPU=$SLURM_NTASKS

ref=%s
fastq=%s
bwa mem -t$CPU $ref -p $fastq | /opt/linux/centos/7.x/x86_64/pkgs/samtools/1.3/bin/samtools view -Sb - | /opt/linux/centos/7.x/x86_64/pkgs/samtools/1.3/bin/samtools sort - -o $fastq\.map-sorted.bam
bamtools filter -in $fastq\.map-sorted.bam -out $fastq\.map-sorted.NM3.bam -tag NM:"<=3"
samtools view $fastq\.map-sorted.NM3.bam | awk '$6!~/S|H/ && $6~/M/' | cut -f3 | uniq |sort | uniq > $fastq\.pacbio_reads.list 
seqtk subseq $ref $fastq\.pacbio_reads.list > $fastq\.pacbio_reads.fa
rm $fastq\.map-sorted.bam $fastq\.map-sorted.NM3.bam $fastq

echo "Done"
''' %(cpu, shell, ref, fq)
    ofile = open(shell, 'w')
    print >> ofile, cmd
    ofile.close()
    job = 'sbatch %s.sh' %(fq)
    print job
    os.system(job)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-r', '--ref')
    parser.add_argument('-c', '--cpu', default='4')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    
    if args.input:
        args.input = os.path.abspath(args.input)   
    if args.ref:
        args.ref   = os.path.abspath(args.ref) 

    for fq in sorted(glob.glob('%s/*.50reads.f*q' %(args.input))):
        write_slurm_shell_bwa(fq, args.ref, args.cpu)        


if __name__ == '__main__':
    main()

