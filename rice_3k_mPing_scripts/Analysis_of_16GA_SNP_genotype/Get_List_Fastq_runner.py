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
        if id_list.has_key(unit[0]):
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
    parser.add_argument('-i', '--input')
    parser.add_argument('--ping')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    ping = os.path.abspath(args.ping)
    home_dir = os.path.split(os.path.realpath(__file__))[0]
    file_list = glob.glob('%s.repeat_*.list' %(args.input))
    cmd = []
    for f in file_list:
        prefix = re.sub(r'.list', r'', f)
        cmd.append('python %s/Get_List_Fastq.py --list %s --fastq %s.te_repeat.ContainingReads.fq' %(home_dir, f, args.input))
        cmd.append('python %s/Split_Fastq2PE.py --input %s' %(home_dir, re.sub(r'.list', r'.fq', f)))
        cmd.append('/rhome/cjinfeng/BigData/software/Velvet/velvet/velveth %s.assembly 31 -shortPaired -fastq -separate %s_1.fq %s_2.fq -short -fastq %s_unpaired.fq' %(prefix, prefix, prefix, prefix))
        cmd.append('/rhome/cjinfeng/BigData/software/Velvet/velvet/velvetg %s.assembly -ins_length 400 -exp_cov 50 -min_contig_lgth 200 -scaffolding yes' %(prefix))
        cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin/bwa mem %s %s.fq > %s.sam' %(ping, prefix, prefix))
        cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools view -bS -o %s.raw.bam %s.sam' %(prefix, prefix)) 
        cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools sort -m 1000000000 %s.raw.bam %s' %(prefix, prefix))
        cmd.append("/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools view -h %s.bam | perl -lane 'if($F[11] =~ /^NM:i:(\d+)$/){print if $1<=2}else{print}'| /opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools view -bS - -o %s.NM2.bam" %(prefix, prefix))
        cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools index %s.bam' %(prefix))
        cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools index %s.NM2.bam' %(prefix))
        cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools mpileup %s.NM2.bam > %s.NM2.mpileup' %(prefix, prefix))

    for c in cmd:
        print c
        os.system(c)

if __name__ == '__main__':
    main()

